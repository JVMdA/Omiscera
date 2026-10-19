"""Multi-omics - MOFA / Factor-based Integration (fifth tile of the Multi-omics submodule).

QUESTION ANSWERED
-----------------
Several omics measured on the SAME samples vary together in a few hidden ways (the treatment, a cell-type mixture, a batch ...). MOFA
(Multi-Omics Factor Analysis; Argelaguet et al. 2018, Mol Syst Biol 14:e8124; MOFA+, Argelaguet et al. 2020, Genome Biol 21:111) finds
a small number of latent FACTORS that explain that variation:
  - the factor VALUES say which samples sit high or low on each factor;
  - the WEIGHTS say which features of each omics carry it;
  - the VARIANCE EXPLAINED (R2) says how much of each omics' variation each factor captures - a factor that explains variance in several
    omics is shared, one that explains it in a single omics is specific to it.
It is descriptive and unsupervised: a factor is a direction of variation, not a pathway and not a cause, and its sign is arbitrary.

WHAT IS HERE AND WHAT IS NOT
----------------------------
The model itself is the authors' own implementation (mofapy2, a pure-Python package, LGPL-3), used as a library and not re-implemented:
variational inference for the Bayesian factor model with automatic relevance determination and sparse (spike-and-slab) weights, Gaussian
likelihoods, missing values handled by the model. What this module adds around it:
  - data preparation (the same drop / log2 / most-variable rules as Correlation-based Integration, so the two tiles read a matrix alike);
  - a DETERMINISTIC factor order (by total variance explained) and sign (the largest weight is positive) - the raw output has neither;
  - STABILITY: the model is refitted with other random seeds and each factor reports how well it reappears (the MOFA authors advise this);
    a factor that does not reappear is flagged, never hidden;
  - tests of each factor against sample groups (Mann-Whitney U for two groups, Kruskal-Wallis for more; Benjamini-Hochberg across factors);
  - enrichment of a factor's WEIGHTS in pathways: a rank-sum test of the member weights against the others (the idea of MOFA's own
    enrichment, after Frost et al. 2015 PCGSE), vectorised here and checked against scipy.

THE VARIANCE-EXPLAINED CONVENTION is mofapy2's own: R2(view, factor) = 1 - sum((Y - z_k w_k')^2) / sum(Y^2) over the centred (and, by default,
scaled-to-equal-variance) data, with the model's prediction set to 0 at missing cells. It is reported as a percentage.

MISSING VALUES AND SAMPLE SIZE
------------------------------
Features missing in more than half of the samples are dropped. MOFA needs enough samples to tell factors apart: fewer than MIN_SAMPLES is
refused and fewer than WARN_SAMPLES gets a warning (the factors, and above all the weights, are unreliable there).
"""

from __future__ import annotations

import contextlib
import io
import math

import numpy as np
import pandas as pd

import Multiomics_correlation as mc

OMICS_ORDER = mc.OMICS_ORDER
OMICS_LABELS = mc.OMICS_LABELS
MIN_SAMPLES = 6
WARN_SAMPLES = 15
MAX_FACTORS = 25
CONVERGENCE_MODES = ("fast", "medium", "slow")
MAX_ITERATIONS = 1000

CONFIG_DEFAULTS = {
    "factors": 10, "convergence": "fast", "top_variance": 2000, "log2": False, "scale_views": True, "max_missing": 0.5,
    "min_r2": 1.0,            # percent: a factor explaining less than this in EVERY view is dropped
    "n_refits": 3, "stability_cut": 0.7, "seed": 42, "group_column": "",
}


# ---------------------------------------------------------------------------
# Preparing the views
# ---------------------------------------------------------------------------

def prepare_views(mats: dict, cfg: dict) -> "tuple[dict, dict]":
    """mats: {omics: samples x features} already aligned to the same samples (Multiomics_correlation.align_samples). Per omics: optional
    log2(x + 1), drop features missing in too many samples or constant, keep the most variable. Returns ({omics: matrix}, {omics: report})."""
    views, reports = {}, {}
    for o in [o for o in OMICS_ORDER if o in mats]:
        out, rep = mc.prepare_matrix(mats[o], bool(cfg["log2"]), float(cfg["max_missing"]), int(cfg["top_variance"]))
        if out.shape[1] < 2:
            raise ValueError(f"{OMICS_LABELS[o]}: fewer than two usable features remain after removing constant features and features "
                             f"missing in more than {float(cfg['max_missing']):.0%} of the samples.")
        views[o], reports[o] = out, rep
    if len(views) < 1:
        raise ValueError("No omics to analyse.")
    return views, reports


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

def _fit_raw(views: dict, cfg: dict, seed: int) -> dict:
    """One mofapy2 fit. Returns the raw arrays: Z (samples x K), W {omics: features x K}, R2 (views x K, fraction), R2 total (views),
    the data the model saw (centred / scaled, per view) and the iteration count."""
    from mofapy2.run.entry_point import entry_point

    names = list(views)
    samples = list(views[names[0]].index)
    data = [[views[o].to_numpy(dtype=float)] for o in names]
    ent = entry_point()
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):      # mofapy2 prints its progress even when 'quiet'
        ent.set_data_options(scale_views=bool(cfg["scale_views"]))
        ent.set_data_matrix(data, likelihoods=["gaussian"] * len(names), views_names=names, groups_names=["all"], samples_names=[samples],
                            features_names=[list(map(str, views[o].columns)) for o in names])
        ent.set_model_options(factors=int(cfg["factors"]), spikeslab_weights=True, ard_factors=True, ard_weights=True)
        ent.set_train_options(iter=MAX_ITERATIONS, convergence_mode=str(cfg["convergence"]), dropR2=float(cfg["min_r2"]) / 100.0,
                              gpu_mode=False, seed=int(seed), verbose=False, quiet=True)
        ent.build()
        ent.run()
    model = ent.model
    exp = model.getExpectations()
    Z = np.asarray(exp["Z"]["E"], dtype=float)
    W = {o: np.asarray(w["E"], dtype=float) for o, w in zip(names, exp["W"])}
    Y = {o: np.asarray(y["E"] if isinstance(y, dict) else y, dtype=float) for o, y in zip(names, exp["Y"])}
    r2 = np.asarray(model.calculate_variance_explained()[0], dtype=float)
    r2_total = np.asarray(model.calculate_variance_explained(total=True)[0], dtype=float)
    n_iter = len(model.train_stats.get("elbo", [])) if hasattr(model, "train_stats") else None
    return {"Z": Z, "W": W, "R2": r2, "R2_total": r2_total, "Y": Y, "samples": samples, "views": names, "n_iterations": n_iter,
            "converged": (n_iter is not None and n_iter < MAX_ITERATIONS)}


def _order_and_sign(raw: dict) -> dict:
    """mofapy2 returns factors in no useful order and with an arbitrary sign. Order by total variance explained (all views), largest first
    (ties by the original index), and flip each factor so that its largest-|weight| feature (over all views) has a positive weight."""
    K = raw["Z"].shape[1]
    total = raw["R2"].sum(axis=0)
    order = sorted(range(K), key=lambda k: (-total[k], k))
    Z = raw["Z"][:, order].copy()
    W = {o: raw["W"][o][:, order].copy() for o in raw["views"]}
    R2 = raw["R2"][:, order].copy()
    for j in range(K):
        best, sign = -1.0, 1.0
        for o in raw["views"]:
            col = W[o][:, j]
            i = int(np.argmax(np.abs(col)))
            if abs(col[i]) > best:
                best, sign = abs(col[i]), (1.0 if col[i] >= 0 else -1.0)
        if sign < 0:
            Z[:, j] *= -1.0
            for o in raw["views"]:
                W[o][:, j] *= -1.0
    out = dict(raw)
    out.update({"Z": Z, "W": W, "R2": R2, "original_index": order})
    return out


def fit_mofa(views: dict, cfg: "dict | None" = None, seed: "int | None" = None) -> dict:
    """Fit MOFA on the prepared views and return the ordered, sign-normalised result."""
    cfg = {**CONFIG_DEFAULTS, **(cfg or {})}
    names = list(views)
    n = views[names[0]].shape[0]
    if n < MIN_SAMPLES:
        raise ValueError(f"MOFA needs at least {MIN_SAMPLES} samples shared by every table; {n} were found.")
    if cfg["convergence"] not in CONVERGENCE_MODES:
        raise ValueError(f"Unknown convergence mode {cfg['convergence']!r}.")
    if not (2 <= int(cfg["factors"]) <= MAX_FACTORS):
        raise ValueError(f"The number of factors must be between 2 and {MAX_FACTORS}.")
    raw = _fit_raw(views, cfg, cfg["seed"] if seed is None else seed)
    if raw["Z"].shape[1] == 0:
        raise ValueError("MOFA found no factor that explains at least "
                         f"{cfg['min_r2']}% of the variance of any omics - the omics do not share (or even each hold) structure the model can use. "
                         "Try fewer features, a lower threshold, or check that the samples are not all alike.")
    return _order_and_sign(raw)


def factor_names(k: int) -> list:
    return [f"Factor{i + 1}" for i in range(k)]


# ---------------------------------------------------------------------------
# Stability across random seeds
# ---------------------------------------------------------------------------

def _abs_corr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """|Pearson r| of every column of a with every column of b (samples x factors each). Constant columns give 0."""
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    na, nb = np.sqrt((a * a).sum(axis=0)), np.sqrt((b * b).sum(axis=0))
    with np.errstate(all="ignore"):
        r = np.abs((a.T @ b) / np.outer(na, nb))
    return np.nan_to_num(r)


def stability(views: dict, cfg: dict, ref: dict) -> dict:
    """Refit with other seeds; each reference factor's stability = the MEDIAN over refits of its best |correlation| with a factor of the refit.
    A factor below cfg['stability_cut'] is flagged unstable. Returns {"per_factor": [...], "n_refits", "cut", "fits_kept": [...]}."""
    cfg = {**CONFIG_DEFAULTS, **(cfg or {})}
    n_refits = int(cfg["n_refits"])
    K = ref["Z"].shape[1]
    best = np.zeros((n_refits, K))
    kept = []
    for r in range(n_refits):
        seed = int(cfg["seed"]) + 1 + r
        try:
            fit = fit_mofa(views, cfg, seed=seed)
        except ValueError:
            # another random start that finds no factor at all is a verdict on the reference factors (none came back), not a reason to lose the run
            kept.append({"seed": seed, "n_factors": 0, "converged": False})
            continue
        kept.append({"seed": seed, "n_factors": int(fit["Z"].shape[1]), "converged": bool(fit["converged"])})
        if fit["Z"].shape[1]:
            best[r] = _abs_corr(ref["Z"], fit["Z"]).max(axis=1)
    med = np.median(best, axis=0) if n_refits else np.full(K, np.nan)
    return {"per_factor": [{"factor": f"Factor{k + 1}", "stability": (None if not n_refits else float(med[k])),
                            "min": (None if not n_refits else float(best[:, k].min())),
                            "unstable": bool(n_refits and med[k] < float(cfg["stability_cut"]))} for k in range(K)],
            "n_refits": n_refits, "cut": float(cfg["stability_cut"]), "fits": kept}


# ---------------------------------------------------------------------------
# Factors against sample groups
# ---------------------------------------------------------------------------

def _bh(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return q
    order = np.argsort(pv)
    adj = np.minimum.accumulate((pv[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def group_tests(Z: np.ndarray, groups: "list | np.ndarray", min_per_group: int = 2) -> dict:
    """Each factor against sample groups: Mann-Whitney U (two groups; effect = rank-biserial correlation, positive when the
    first-listed group is higher) or Kruskal-Wallis (more; effect = eta-squared based on H, (H - k + 1) / (n - k), Tomczak & Tomczak 2014), Benjamini-Hochberg across the factors.
    Samples without a group label ('' / None) are left out. Groups with fewer than min_per_group samples are dropped."""
    from scipy import stats

    g = np.array([("" if x is None else str(x)) for x in groups], dtype=object)
    labels = [x for x in dict.fromkeys(g.tolist()) if x != "" and (g == x).sum() >= min_per_group]
    out = {"groups": labels, "counts": {x: int((g == x).sum()) for x in labels}, "rows": [], "test": None}
    if len(labels) < 2:
        out["note"] = "Fewer than two groups with at least two samples each, so no group test was run."
        return out
    keep = np.isin(g, labels)
    gk = g[keep]
    out["test"] = "Mann-Whitney U" if len(labels) == 2 else "Kruskal-Wallis"
    pvals, eff = [], []
    for k in range(Z.shape[1]):
        z = Z[keep, k]
        samples = [z[gk == x] for x in labels]
        if len(labels) == 2:
            res = stats.mannwhitneyu(samples[0], samples[1], alternative="two-sided")
            n1, n2 = len(samples[0]), len(samples[1])
            eff.append(float(2.0 * res.statistic / (n1 * n2) - 1.0))
        else:
            res = stats.kruskal(*samples)
            n = int(keep.sum())
            eff.append(float((res.statistic - len(labels) + 1) / (n - len(labels))) if n > len(labels) else float("nan"))
        pvals.append(float(res.pvalue))
    q = _bh(pvals)
    for k in range(Z.shape[1]):
        row = {"Factor": f"Factor{k + 1}", "Test": out["test"], "p-value": pvals[k], "FDR": float(q[k]), "Effect": eff[k]}
        for x in labels:
            row[f"Median ({x})"] = float(np.median(Z[keep, k][gk == x]))
        out["rows"].append(row)
    return out


# ---------------------------------------------------------------------------
# Enrichment of a factor's weights in pathways
# ---------------------------------------------------------------------------

def weight_set_enrichment(weights: np.ndarray, keys: list, term_map: dict, min_size: int = 10, max_size: int = 500, mode: str = "signed") -> pd.DataFrame:
    """Are a pathway's members high (or low) among the weights of a factor?
    weights: one value per feature; keys: the identifier each feature has in the pathway database (None / '' for a feature that has none -
    it is left out of the test). term_map: {term id: {"name", "genes": set of keys}}.
    mode 'signed': two-sided rank-sum test on the signed weights (members vs the other features); AUC > 0.5 = the members' weights are higher;
    mode 'absolute': one-sided test on |weights| (is the pathway involved in the factor, whatever the sign).
    The test is Mann-Whitney U with the normal approximation, tie-corrected and with continuity correction (the same as scipy's
    asymptotic Mann-Whitney), Benjamini-Hochberg over the pathways tested. Returns a DataFrame sorted by FDR."""
    from scipy import stats

    if mode not in ("signed", "absolute"):
        raise ValueError("mode must be 'signed' or 'absolute'.")
    w = np.asarray(weights, dtype=float)
    ok = np.array([bool(k) for k in keys]) & np.isfinite(w)
    idx_of: dict = {}
    for i in np.nonzero(ok)[0]:
        idx_of.setdefault(keys[i], []).append(int(i))
    cols = ["Term ID", "Term Name", "Members", "AUC", "z", "p-value", "FDR", "Direction", "Top members"]
    sub = np.nonzero(ok)[0]
    if sub.size < 2 * min_size:
        return pd.DataFrame(columns=cols)
    vals = np.abs(w[sub]) if mode == "absolute" else w[sub]
    ranks = stats.rankdata(vals, method="average")
    pos = {int(i): j for j, i in enumerate(sub)}          # feature index -> position in the tested vector
    N = len(sub)
    _, counts = np.unique(vals, return_counts=True)
    tie = float((counts ** 3 - counts).sum())
    rows = []
    for tid, info in term_map.items():
        mem = sorted({pos[i] for k in info["genes"] for i in idx_of.get(k, ())})
        n1 = len(mem)
        n2 = N - n1
        if n1 < min_size or n1 > max_size or n2 < min_size:
            continue
        u1 = float(ranks[mem].sum() - n1 * (n1 + 1) / 2.0)
        mu = n1 * n2 / 2.0
        sigma = math.sqrt(n1 * n2 / 12.0 * ((N + 1) - tie / (N * (N - 1))))
        if sigma <= 0:
            continue
        if mode == "signed":
            z = (abs(u1 - mu) - 0.5) / sigma if u1 != mu else 0.0
            z = math.copysign(max(z, 0.0), u1 - mu)
            p = float(2.0 * stats.norm.sf(abs(z)))
        else:
            z = (u1 - mu - 0.5) / sigma
            p = float(stats.norm.sf(z))
        top = sorted(mem, key=lambda j: -abs(vals[j]))[:5]
        rows.append({"Term ID": tid, "Term Name": info["name"], "Members": n1, "AUC": u1 / (n1 * n2), "z": float(z), "p-value": min(max(p, 0.0), 1.0),
                     "Direction": ("higher weights" if u1 > mu else "lower weights") if mode == "signed" else "larger |weight|",
                     "Top members": ", ".join(str(keys[int(sub[j])]) for j in top)})
    df = pd.DataFrame(rows, columns=[c for c in cols if c != "FDR"])
    if df.empty:
        return pd.DataFrame(columns=cols)
    df["FDR"] = _bh(df["p-value"].to_numpy())
    return df.sort_values(by=["FDR", "p-value", "Term ID"], kind="stable").reset_index(drop=True)[cols]


def collapse_by_key(weights: np.ndarray, keys: list) -> "tuple[np.ndarray, list]":
    """Several features can stand for one gene (probes, isoforms, protein groups). For the pathway test a gene must count once: keep, for each
    key, the feature with the largest |weight| (its sign kept). Features without a key are left as they are. Order = first appearance."""
    w = np.asarray(weights, dtype=float)
    best: dict = {}
    for i, k in enumerate(keys):
        if k and (k not in best or abs(w[i]) > abs(w[best[k]])):
            best[k] = i
    keep = [i for i, k in enumerate(keys) if (not k) or best[k] == i]
    return w[keep], [keys[i] for i in keep]


# ---------------------------------------------------------------------------
# Putting the whole analysis together
# ---------------------------------------------------------------------------

def weights_table(result: dict, views: dict, labels: dict) -> dict:
    """{omics: DataFrame(feature id, label, Factor1.. weights)}."""
    K = result["Z"].shape[1]
    out = {}
    for o in result["views"]:
        df = pd.DataFrame(result["W"][o], columns=factor_names(K))
        df.insert(0, "Label", [labels.get(o, {}).get(str(c), str(c)) for c in views[o].columns])
        df.insert(0, "Feature", [str(c) for c in views[o].columns])
        out[o] = df
    return out


def run_mofa(mats: dict, cfg: "dict | None" = None, groups: "dict | None" = None, labels: "dict | None" = None, progress=None) -> dict:
    """The whole analysis on aligned matrices. groups: {sample: group label} for the interpretation only (the model is fitted on all samples
    together). labels: {omics: {feature id: display label (gene symbol ...)}}.
    Returns {"config", "samples", "omics", "views": {omics: matrix used}, "reports", "fit", "factors", "r2", "r2_total", "stability",
    "weights", "group_tests", "warnings"}."""
    cfg = {**CONFIG_DEFAULTS, **(cfg or {})}
    labels = labels or {}

    def _p(f, text):
        if progress:
            try:
                progress(f, text)
            except Exception:
                pass

    n = next(iter(mats.values())).shape[0]
    if n < MIN_SAMPLES:
        raise ValueError(f"MOFA needs at least {MIN_SAMPLES} samples shared by every table; {n} were found.")
    warnings = []
    if n < WARN_SAMPLES:
        warnings.append(f"Only {n} samples: with fewer than {WARN_SAMPLES} MOFA's factors, and above all its weights, are unreliable. Treat the "
                        "result as exploratory and look at the stability column.")
    _p(0.05, "Preparing the data...")
    views, reports = prepare_views(mats, cfg)
    if len(views) < 2:
        warnings.append("Only one omics remained, so there is nothing shared to find: MOFA then reduces to a factor analysis of that omics alone.")
    _p(0.15, "Fitting the MOFA model...")
    fit = fit_mofa(views, cfg)
    K = fit["Z"].shape[1]
    if not fit["converged"]:
        warnings.append(f"The model did not converge within {MAX_ITERATIONS} iterations; use 'medium' or 'slow' convergence for a more careful fit.")
    _p(0.55, "Refitting with other seeds to check stability...")
    stab = stability(views, cfg, fit) if int(cfg["n_refits"]) > 0 else {"per_factor": [{"factor": f"Factor{k + 1}", "stability": None, "min": None, "unstable": False} for k in range(K)],
                                                                         "n_refits": 0, "cut": float(cfg["stability_cut"]), "fits": []}
    unstable = [p["factor"] for p in stab["per_factor"] if p["unstable"]]
    if unstable:
        warnings.append(f"{', '.join(unstable)} did not reappear (|r| below {stab['cut']}) when the model was refitted with other seeds: "
                        "read it as unstable, not as a finding.")
    _p(0.9, "Testing the factors against the sample groups...")
    gt = group_tests(fit["Z"], [(groups or {}).get(s, "") for s in fit["samples"]]) if groups else {"groups": [], "counts": {}, "rows": [], "test": None,
                                                                                                    "note": "No sample groups were given."}
    names = factor_names(K)
    r2 = pd.DataFrame(100.0 * fit["R2"], index=[OMICS_LABELS[o] for o in fit["views"]], columns=names)
    r2_total = {o: float(100.0 * fit["R2_total"][i]) for i, o in enumerate(fit["views"])}
    _p(1.0, "Done")
    return {"config": cfg, "samples": fit["samples"], "omics": list(fit["views"]), "views": views, "reports": reports, "fit": fit,
            "factors": pd.DataFrame(fit["Z"], index=fit["samples"], columns=names), "r2": r2, "r2_total": r2_total, "stability": stab,
            "weights": weights_table(fit, views, labels), "group_tests": gt, "warnings": warnings}
