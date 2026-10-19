"""Multi-omics - Correlation-based Integration (third tile of the Multi-omics submodule).

QUESTIONS ANSWERED
------------------
Given a SAMPLE x FEATURE quantitative table per omics - the same samples in every file - two questions:

  1. SAME-GENE CONCORDANCE. For every gene present in both the transcriptomics and the proteomics table: how well
     does its transcript level track its protein level across the samples (rho, p, FDR)? The classic RNA-versus-
     protein question (Vogel & Marcotte 2012, Nat Rev Genet 13:227; Liu, Beyer & Aebersold 2016, Cell 165:535).
  2. CROSS-OMICS FEATURE PAIRS. Every feature of one omics against every feature of another (RNA-protein,
     RNA-metabolite, protein-metabolite): the strongest pairs after FDR, as a table, a heatmap and a network.

It is descriptive: a correlation across samples is an association, not a mechanism, and where the samples fall into
groups (control / treated) a shared group effect makes unrelated features correlate. The README says so.

THE STATISTICS
--------------
  - Correlation: Spearman (default; Pearson on average ranks, each PAIR ranked over the samples both features have - the
    definition R and pandas use with pairwise-complete data) or Pearson. Pairwise-complete samples, each pair with its
    own n; a pair needs at least MIN_SAMPLES.
  - p-value: t = r * sqrt((n - 2) / (1 - r^2)), two-sided, n - 2 degrees of freedom (the usual test of r = 0; for
    Spearman it is the standard large-sample approximation, so it is unreliable at tiny n - fewer than WARN_SAMPLES
    shared samples is warned about, fewer than MIN_SAMPLES refused).
  - Multiple testing: Benjamini-Hochberg, one family per block (RNA-protein, RNA-metabolite, protein-metabolite are three
    families; the same-gene analysis is its own). The number of tests per block and the number of reported pairs that
    are expected to be false (the sum of their q-values) are reported, not hidden.

THIS ENVIRONMENT
----------------
numpy.corrcoef / numpy.cov and scipy.stats.spearmanr HARD-CRASH the interpreter here (a BLAS fault), so nothing in this
module uses them or any matrix product: every correlation is an elementwise multiply-and-sum, one feature of A against
all features of B at a time.

DESIGN RULES (each one a decision, tested in tests/test_multiomics_correlation.py)
----------------------------------------------------------------------------------
1. Samples are matched by their id (whitespace trimmed, case optional) and the analysis uses ONLY the samples every
   omics in the run has. Which were dropped, per omics, is reported - never silently. Fewer than MIN_SAMPLES shared
   samples is refused, with a hint when the other orientation of a file would match.
2. A feature with too many missing values (> MAX_MISSING_FRACTION), or constant, cannot be correlated and is dropped and
   counted. Values are taken as already normalised; the optional log2(x + 1) is for raw counts / intensities.
3. The cross-omics blocks keep the top-N features by variance per omics (a pair count of n^2 has to stay tractable);
   the same-gene analysis uses ALL genes (it is only one test per gene).
4. Several features of one omics that resolve to the same gene are collapsed by their mean before the same-gene test.
5. A pair is REPORTED when its q-value is at most the FDR cutoff AND |rho| is at least the rho cutoff. The BH family is
   every tested pair of the block, never only the reported ones.
"""

from __future__ import annotations

import io
import math

import numpy as np
import pandas as pd

MIN_SAMPLES = 5
WARN_SAMPLES = 10
MAX_MISSING_FRACTION = 0.2
METHODS = ("spearman", "pearson")
METHOD_DEFAULT = "spearman"
OMICS_ORDER = ("transcriptomics", "proteomics", "metabolomics")
OMICS_LABELS = {"transcriptomics": "Transcriptomics", "proteomics": "Proteomics", "metabolomics": "Metabolomics"}
MAX_PAIRS_PER_BLOCK = 16_000_000
MAX_REPORTED_ROWS = 20_000

CONFIG_DEFAULTS = {
    "method": METHOD_DEFAULT, "log2": False, "top_variance": 2000, "fdr": 0.05, "rho_cut": 0.6,
    "min_samples": MIN_SAMPLES, "case_insensitive": False, "heat_k": 40, "net_k": 60,
}

# A text column named with one of these words labels SAMPLES (MetaboAnalyst's sample-in-rows format calls it "Label").
GROUP_COLUMN_WORDS = {"group", "groups", "condition", "conditions", "treatment", "phenotype", "cohort", "genotype", "arm",
                      "label", "labels"}
# These also name sample groups but just as often annotate FEATURES (a compound's class, a gene's type): they count only
# as "sample class" / "sample type" style names.
GROUP_COLUMN_WEAK_WORDS = {"class", "type", "status"}
# A column named with any of these annotates features, never samples (gene_biotype, Compound class, Protein names ...).
FEATURE_COLUMN_WORDS = {"gene", "genes", "protein", "proteins", "compound", "compounds", "metabolite", "metabolites",
                        "feature", "features", "transcript", "transcripts", "biotype", "chemical", "pathway", "pathways",
                        "superclass", "subclass", "kegg", "hmdb", "chebi", "description", "annotation", "ontology",
                        "symbol", "accession"}


def _group_like_column(name: str, cells: pd.Series) -> bool:
    """Does this set-aside text column look like a per-sample group label (so each ROW is a sample)?"""
    import re
    tokens = {t for t in re.split(r"[^a-z0-9]+", str(name).lower()) if t}
    if not tokens or tokens & FEATURE_COLUMN_WORDS:
        return False
    if not (tokens & GROUP_COLUMN_WORDS or (tokens & GROUP_COLUMN_WEAK_WORDS and "sample" in tokens)):
        return False
    vals = cells[(cells != "") & (cells.str.lower() != "nan")]
    return len(vals) >= 2 and vals.nunique() < len(vals)      # a group labels several samples: some value repeats


# ---------------------------------------------------------------------------
# Reading a matrix upload
# ---------------------------------------------------------------------------

def _decode(raw: bytes) -> str:
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def _read_table(raw: bytes) -> pd.DataFrame:
    text = _decode(raw)
    if not text.strip():
        raise ValueError("The file is empty.")
    for sep in (None, ",", "\t", ";"):
        try:
            kw = {"sep": sep, "engine": "python"} if sep is None else {"sep": sep}
            df = pd.read_csv(io.StringIO(text), dtype=str, **kw)
        except Exception:
            continue
        if df.shape[1] >= 3 and df.shape[0] >= 2:
            return df
    raise ValueError("Could not read the file as a table with an identifier column and at least two data columns "
                     "(CSV / TSV expected).")


def _uniquify(labels) -> "tuple[list, int]":
    seen: dict = {}
    out, n_dup = [], 0
    for lab in labels:
        base = str(lab)
        if base in seen:
            seen[base] += 1
            n_dup += 1
            out.append(f"{base} ({seen[base]})")
        else:
            seen[base] = 1
            out.append(base)
    return out, n_dup


def parse_matrix(raw_bytes: bytes, orientation: str = "auto") -> dict:
    """Read a quantitative table. The first column holds identifiers; every other column that is mostly numbers is data,
    the rest (a Group column, say) is set aside and counted.

    orientation: "columns" (samples in the columns, features in the rows - the usual expression-matrix layout), "rows"
    (samples in the rows), or "auto": a group-like text column means samples are in the rows (one Group label per sample);
    otherwise the axis with fewer labels is taken to be the samples (an RNA matrix has thousands of rows, a handful of
    sample columns).

    Returns {"df": samples x features (float), "meta": {...}}. meta["alt_samples"] is the label set of the OTHER axis, so a
    sample mismatch can suggest the file is the wrong way round."""
    if orientation not in ("auto", "columns", "rows"):
        raise ValueError("orientation must be auto, columns or rows.")
    table = _read_table(raw_bytes)
    table.columns = [str(c).strip() for c in table.columns]
    id_name = table.columns[0]
    ids = table.iloc[:, 0].astype(str).str.strip()
    keep_rows = ids != ""
    table, ids = table.loc[keep_rows], ids[keep_rows]
    data, set_aside, group_cols = {}, [], []
    for col in table.columns[1:]:
        cells = table[col].astype(str).str.strip()
        nonblank = cells[(cells != "") & (cells.str.lower() != "nan")]
        numeric = pd.to_numeric(cells, errors="coerce")
        if len(nonblank) and numeric.notna().sum() >= 0.5 * len(nonblank):
            data[col] = numeric.to_numpy(dtype=float)
        else:
            set_aside.append(col)
            if _group_like_column(col, cells):
                group_cols.append(col)
    if len(data) < 2:
        raise ValueError("Fewer than two numeric data columns were found.")
    raw_df = pd.DataFrame(data, index=ids.to_numpy())
    reason = "set by you"
    chosen = orientation
    if orientation == "auto":
        if group_cols:
            chosen, reason = "rows", f"a group-like column ('{group_cols[0]}') labels each row, so each row is a sample"
        elif raw_df.shape[0] > raw_df.shape[1]:
            chosen, reason = "columns", (f"{raw_df.shape[0]} rows against {raw_df.shape[1]} data columns: the columns "
                                         "are the samples")
        else:
            chosen, reason = "rows", (f"{raw_df.shape[0]} rows against {raw_df.shape[1]} data columns: the rows "
                                      "are the samples")
    if chosen == "columns":
        samples = list(raw_df.columns)
        feat_labels = list(raw_df.index)
        values = raw_df.to_numpy().T
        alt = [str(x) for x in raw_df.index]
    else:
        samples = list(raw_df.index)
        feat_labels = list(raw_df.columns)
        values = raw_df.to_numpy()
        alt = [str(x) for x in raw_df.columns]
    samples = [str(s).strip() for s in samples]
    if len(set(samples)) != len(samples):
        dup = sorted({s for s in samples if samples.count(s) > 1})[:5]
        raise ValueError(f"Duplicate sample ids ({', '.join(dup)}). Each sample must appear once - if the samples are in "
                         "the other direction, set the orientation.")
    feat_ids, n_dup = _uniquify(feat_labels)
    df = pd.DataFrame(values, index=samples, columns=feat_ids)
    n_before = df.shape
    df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")
    if df.shape[0] < 3 or df.shape[1] < 2:
        raise ValueError("After removing empty rows and columns fewer than 3 samples or 2 features remain.")
    # A set-aside text column with a few repeated values, in a table whose rows are the SAMPLES, is a per-sample label (Group, Condition ...).
    # It plays no part in the correlation but MOFA uses it to colour and test its factors, so it is kept: {column: {sample: value}}.
    group_columns: dict = {}
    if chosen == "rows":
        for col in set_aside:
            vals = table[col].astype(str).str.strip().to_numpy()
            mp = {str(s).strip(): v for s, v in zip(ids.to_numpy(), vals) if v not in ("", "nan")}
            distinct = set(mp.values())
            if 2 <= len(distinct) <= max(2, len(mp) // 2):
                group_columns[col] = {s: mp[s] for s in df.index if s in mp}
    meta = {"orientation": chosen, "orientation_reason": reason, "id_column": id_name, "group_columns": group_columns,
            "n_samples": int(df.shape[0]), "n_features": int(df.shape[1]), "set_aside_columns": set_aside,
            "n_repeated_feature_ids": n_dup, "n_empty_dropped": int(n_before[0] - df.shape[0] + n_before[1] - df.shape[1]),
            "first_samples": list(df.index[:5]), "alt_samples": alt[:5000],
            # a repeated id was renamed "X (2)": its identifier for resolution is still "X"
            "original_ids": {u: str(o) for u, o in zip(feat_ids, feat_labels) if u != str(o)}}
    return {"df": df, "meta": meta}


# ---------------------------------------------------------------------------
# Cleaning, sample matching
# ---------------------------------------------------------------------------

def prepare_matrix(df: pd.DataFrame, log2: bool = False, max_missing: float = MAX_MISSING_FRACTION,
                   top_variance: "int | None" = None) -> "tuple[pd.DataFrame, dict]":
    """Drop the features that cannot be correlated (too many missing, constant), optionally log2(x + 1), and optionally
    keep the top-N by variance. Returns (matrix, report)."""
    report: dict = {"n_in": int(df.shape[1])}
    out = df.astype(float)
    if log2:
        if (out.to_numpy() < 0).any():
            raise ValueError("log2(x + 1) needs non-negative values; this table has negative ones (already log-scaled?).")
        out = np.log2(out + 1.0)
    miss = out.isna().mean(axis=0)
    too_missing = miss > max_missing
    report["n_too_many_missing"] = int(too_missing.sum())
    out = out.loc[:, ~too_missing.to_numpy()]
    arr = out.to_numpy()
    n_obs = np.sum(~np.isnan(arr), axis=0)
    with np.errstate(all="ignore"):
        mean = np.nanmean(arr, axis=0)
        var = np.nanmean((arr - mean) ** 2, axis=0)
    constant = (n_obs < 3) | ~(var > 1e-12 * np.maximum(1.0, np.abs(mean) ** 2)) | ~np.isfinite(var)
    report["n_constant"] = int(constant.sum())
    out = out.loc[:, ~constant]
    var = var[~constant]
    report["n_kept_before_variance"] = int(out.shape[1])
    if top_variance and out.shape[1] > top_variance:
        order = np.argsort(-var, kind="mergesort")[:top_variance]
        out = out.iloc[:, np.sort(order)]
    report["n_out"] = int(out.shape[1])
    report["n_variance_filtered"] = report["n_kept_before_variance"] - report["n_out"]
    return out, report


def align_samples(mats: dict, case_insensitive: bool = False, min_samples: int = MIN_SAMPLES) -> dict:
    """Keep the samples every matrix has. mats: {omics: samples x features}. Returns {"samples", "mats", "report",
    "warnings"}; raises ValueError (naming a likelier cause) below min_samples."""
    def norm(s):
        s = str(s).strip()
        return s.lower() if case_insensitive else s

    keyed: dict = {}
    for o, df in mats.items():
        keys = [norm(s) for s in df.index]
        if len(set(keys)) != len(keys):
            raise ValueError(f"{OMICS_LABELS.get(o, o)}: two samples have the same id once whitespace"
                             f"{' and case' if case_insensitive else ''} is ignored.")
        keyed[o] = dict(zip(keys, df.index))
    shared_keys = None
    for o, k in keyed.items():
        shared_keys = set(k) if shared_keys is None else shared_keys & set(k)
    first = next(iter(mats))
    ordered = [k for k in keyed[first] if k in shared_keys]
    report, warnings = {}, []
    for o, df in mats.items():
        dropped = [str(s) for s in df.index if norm(s) not in shared_keys]
        report[o] = {"n_samples": int(df.shape[0]), "n_used": len(ordered), "dropped": dropped}
    if len(ordered) < min_samples:
        hint = ""
        for o, df in mats.items():
            alt = getattr(df, "attrs", {}).get("alt_samples")
            if alt:
                others = [set(norm(s) for s in m.index) for oo, m in mats.items() if oo != o]
                if others and len(set(norm(s) for s in alt).intersection(*others)) >= min_samples:
                    hint = (f" The other orientation of the {OMICS_LABELS.get(o, o)} file would share "
                            f"{len(set(norm(s) for s in alt).intersection(*others))} samples with the others: "
                            "its samples may be in the other direction - change its orientation.")
                    break
        raise ValueError(f"Only {len(ordered)} sample(s) are in every file; at least {min_samples} are needed." + hint)
    if len(ordered) < WARN_SAMPLES:
        warnings.append(f"Only {len(ordered)} shared samples: the p-values rest on an approximation that is unreliable at "
                        f"this size - treat the results as exploratory.")
    out = {o: df.loc[[keyed[o][k] for k in ordered]].copy() for o, df in mats.items()}
    names = [str(keyed[first][k]).strip() for k in ordered]
    for o, df in out.items():
        df.index = names
    return {"samples": names, "mats": out, "report": report, "warnings": warnings}


# ---------------------------------------------------------------------------
# Correlation
# ---------------------------------------------------------------------------

def _rankdata(x: np.ndarray) -> np.ndarray:
    """Average ranks (ties share the mean of their ranks), 1-based."""
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt)
    return ((csum - cnt + 1) + csum)[inv] / 2.0


def _rank_columns(M: np.ndarray) -> np.ndarray:
    """Average ranks of every column over the values it HAS (NaN stays NaN), vectorised over the columns: sort each
    column (NaN sorted last as +inf), give each run of equal values the mean of its positions, put the ranks back."""
    X = np.asarray(M, dtype=float)
    if X.size == 0:
        return np.full(X.shape, np.nan)
    n, m = X.shape
    nan = np.isnan(X)
    Xf = np.where(nan, np.inf, X)
    order = np.argsort(Xf, axis=0, kind="mergesort")
    xs = np.take_along_axis(Xf, order, axis=0)
    pos = np.arange(1, n + 1, dtype=float)[:, None]
    start = np.ones((n, m), dtype=bool)
    start[1:] = xs[1:] != xs[:-1]
    end = np.ones((n, m), dtype=bool)
    end[:-1] = xs[1:] != xs[:-1]
    first = np.maximum.accumulate(np.where(start, pos, 0.0), axis=0)
    last = np.minimum.accumulate(np.where(end, pos, np.inf)[::-1], axis=0)[::-1]
    out = np.empty((n, m))
    np.put_along_axis(out, order, (first + last) / 2.0, axis=0)
    out[nan] = np.nan
    return out


def _pearson_sums(ra: np.ndarray, rb: np.ndarray, valid: np.ndarray):
    """Pearson r of column pairs over the rows marked valid. ra, rb, valid: n x q. Returns (r, n)."""
    a0, b0 = np.where(valid, ra, 0.0), np.where(valid, rb, 0.0)
    npair = valid.sum(0)
    sx, sy = a0.sum(0), b0.sum(0)
    with np.errstate(all="ignore"):
        num = npair * (a0 * b0).sum(0) - sx * sy
        den = np.sqrt((npair * (a0 * a0).sum(0) - sx * sx) * (npair * (b0 * b0).sum(0) - sy * sy))
        r = num / den
    r[~(den > 0)] = np.nan
    return r, npair


def _spearman_pairwise(A: np.ndarray, B: np.ndarray, progress=None):
    """Exact Spearman for every A column x B column when values are missing: each pair is ranked over ITS OWN shared
    samples (as R's cor(use = "pairwise.complete.obs", method = "spearman") and pandas do), never over the samples one
    feature has alone - that shortcut moves rho by up to ~0.1 at 12 samples with one missing value.
    For a feature a, all of B at once:
      - B's ranks over (B's samples AND a's samples): one vectorised column ranking per distinct missing pattern of a;
      - a's ranks over (a's samples AND b's samples) for every b: rank_v(a_s) = sum over k in v of w(k, s) + 1/2, with
        w(k, s) = [a_k < a_s] + [a_k = a_s] / 2 - an einsum over the samples (no BLAS matrix product, which crashes here)."""
    n, p = A.shape
    q = B.shape[1]
    MA, MB = ~np.isnan(A), ~np.isnan(B)
    MBf = MB.astype(float)
    b_cols_missing = np.nonzero(~MB.all(0))[0]
    R = np.full((p, q), np.nan, dtype=np.float64)
    N = np.zeros((p, q), dtype=np.int16)
    RB_own = _rank_columns(B)
    by_mask = sorted(range(p), key=lambda i: MA[:, i].tobytes())     # same pattern -> consecutive -> one B ranking
    last_key, RB = None, None
    for done, i in enumerate(by_mask):
        if progress and done % 100 == 0:
            progress(done / max(p, 1))
        a, ma = A[:, i], MA[:, i]
        key = ma.tobytes()
        if key != last_key:
            RB = RB_own if ma.all() else _rank_columns(np.where(ma[:, None], B, np.nan))
            last_key = key
        ra_own = _rank_columns(a[:, None])[:, 0]
        RA = np.repeat(ra_own[:, None], q, axis=1)
        if b_cols_missing.size:
            a0 = np.where(ma, a, np.nan)
            W = (a0[:, None] < a0[None, :]) + 0.5 * (a0[:, None] == a0[None, :])     # NaN compares False: no weight
            RA[:, b_cols_missing] = np.einsum("ks,kj->sj", W, MBf[:, b_cols_missing], optimize=False) + 0.5
        r, npair = _pearson_sums(RA, RB, ma[:, None] & MB)
        R[i] = r
        N[i] = npair
    return R, N


def _check_method(method: str) -> None:
    if method not in METHODS:
        raise ValueError(f"Unknown correlation method {method!r}.")


def correlate_all(A: np.ndarray, B: np.ndarray, method: str = METHOD_DEFAULT, min_n: int = MIN_SAMPLES, progress=None):
    """Every feature of A (columns) against every feature of B. Returns (R, N): the correlation (float64 - a float32 r loses the accuracy of the p-value of a very strong correlation; NaN where the
    pair has fewer than min_n common samples or no spread) and the number of samples each pair used (int16).
    Elementwise numpy only - never a matrix product."""
    _check_method(method)
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    if method == "spearman":
        if np.isnan(A).any() or np.isnan(B).any():
            R, N = _spearman_pairwise(A, B, progress)
            np.clip(R, -1.0, 1.0, out=R)
            R[N < min_n] = np.nan
            return R, N
        A, B = _rank_columns(A), _rank_columns(B)
    n, p = A.shape
    q = B.shape[1]
    R = np.full((p, q), np.nan, dtype=np.float64)
    N = np.zeros((p, q), dtype=np.int16)
    Bm = ~np.isnan(B)
    B0 = np.where(Bm, B, 0.0)
    B2 = B0 * B0
    full = (not np.isnan(A).any()) and bool(Bm.all())
    if full:
        with np.errstate(all="ignore"):
            za = (A - A.mean(0)) / A.std(0, ddof=1)
            zb = (B - B.mean(0)) / B.std(0, ddof=1)
        for i in range(p):
            if progress and i % 100 == 0:
                progress(i / max(p, 1))
            R[i] = (za[:, i][:, None] * zb).sum(0) / (n - 1)
        N[:] = n
        R[~np.isfinite(R)] = np.nan
    else:
        for i in range(p):
            if progress and i % 100 == 0:
                progress(i / max(p, 1))
            a = A[:, i]
            ma = ~np.isnan(a)
            a0 = np.where(ma, a, 0.0)
            valid = ma[:, None] & Bm
            npair = valid.sum(0)
            sx = (a0[:, None] * valid).sum(0)
            sy = (B0 * ma[:, None]).sum(0)
            sxx = ((a0 * a0)[:, None] * valid).sum(0)
            syy = (B2 * ma[:, None]).sum(0)
            sxy = (a0[:, None] * B0).sum(0)
            with np.errstate(all="ignore"):
                num = npair * sxy - sx * sy
                den = np.sqrt((npair * sxx - sx * sx) * (npair * syy - sy * sy))
                r = num / den
            r[~(den > 0)] = np.nan
            R[i] = r
            N[i] = npair
    np.clip(R, -1.0, 1.0, out=R)
    R[N < min_n] = np.nan
    return R, N


def paired_correlation(A: np.ndarray, B: np.ndarray, method: str = METHOD_DEFAULT, min_n: int = MIN_SAMPLES):
    """Column j of A against column j of B (the same gene's transcript and protein). Returns (r, n) vectors."""
    _check_method(method)
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    if A.shape != B.shape:
        raise ValueError("The two matrices must have the same shape.")
    valid = ~np.isnan(A) & ~np.isnan(B)
    if method == "spearman":
        # each pair ranked over the samples BOTH have (see _spearman_pairwise)
        A = _rank_columns(np.where(valid, A, np.nan))
        B = _rank_columns(np.where(valid, B, np.nan))
    r, n = _pearson_sums(A, B, valid)
    r = np.clip(r, -1.0, 1.0)
    r[n < min_n] = np.nan
    return r, n


def pvalues_from_r(R, N) -> np.ndarray:
    """Two-sided p of r = 0: t = r * sqrt((n - 2) / (1 - r^2)) on n - 2 degrees of freedom. NaN where r is NaN or n < 3."""
    from scipy.stats import t as _t

    r = np.asarray(R, dtype=np.float64)
    n = np.asarray(N, dtype=np.float64)
    out = np.full(r.shape, np.nan)
    ok = np.isfinite(r) & (n >= 3)
    rr, nn = r[ok], n[ok]
    with np.errstate(all="ignore"):
        tt = rr * np.sqrt((nn - 2.0) / (1.0 - rr * rr))
    p = 2.0 * _t.sf(np.abs(tt), nn - 2.0)
    p[np.abs(rr) >= 1.0] = 0.0
    out[ok] = p
    return out


def bh_qvalues(p_all, p_query) -> np.ndarray:
    """Benjamini-Hochberg q-values of p_query inside the family p_all (finite values only). Equivalent to the usual BH
    over p_all, evaluated for the queried p-values without ranking every pair a second time."""
    p = np.asarray(p_all, dtype=np.float64).ravel()
    p = p[np.isfinite(p)]
    q_in = np.asarray(p_query, dtype=np.float64)
    out = np.full(q_in.shape, np.nan)
    m = p.size
    if m == 0:
        return out
    s = np.sort(p)
    vals = s * m / np.arange(1, m + 1)
    suffix = np.minimum(np.minimum.accumulate(vals[::-1])[::-1], 1.0)
    ok = np.isfinite(q_in)
    idx = np.clip(np.searchsorted(s, q_in[ok], side="right") - 1, 0, m - 1)
    out[ok] = suffix[idx]
    return out


def bh_threshold(p_all, fdr: float) -> "tuple[float, int]":
    """(largest p that is significant at this FDR under BH, how many p-values are). (0.0, 0) when none."""
    p = np.asarray(p_all, dtype=np.float64).ravel()
    p = p[np.isfinite(p)]
    m = p.size
    if m == 0:
        return 0.0, 0
    s = np.sort(p)
    ok = np.nonzero(s <= np.arange(1, m + 1) * fdr / m)[0]
    if ok.size == 0:
        return 0.0, 0
    k = int(ok[-1]) + 1
    return float(s[k - 1]), k


# ---------------------------------------------------------------------------
# Same-gene concordance
# ---------------------------------------------------------------------------

def collapse_by_gene(df: pd.DataFrame, gene_of: dict) -> "tuple[pd.DataFrame, dict]":
    """Samples x genes: features that resolve to the same gene are averaged (mean over the features that have a value);
    features with no gene are left out. gene_of: {feature id: GENE (upper-case) or ''}.
    Returns (matrix, info) with info["members"] = {gene: [feature ids]}."""
    members: dict = {}
    n_unmapped = 0
    for col in df.columns:
        g = str(gene_of.get(col, "") or "").strip().upper()
        if not g:
            n_unmapped += 1
            continue
        members.setdefault(g, []).append(col)
    cols, names = [], []
    for g in sorted(members):
        feats = members[g]
        cols.append(df[feats[0]].to_numpy() if len(feats) == 1 else np.nanmean(df[feats].to_numpy(), axis=1))
        names.append(g)
    out = pd.DataFrame(np.column_stack(cols) if cols else np.empty((df.shape[0], 0)), index=df.index, columns=names)
    return out, {"members": members, "n_in": int(df.shape[1]), "n_genes": len(names), "n_unmapped": n_unmapped,
                 "n_collapsed_genes": sum(1 for v in members.values() if len(v) > 1)}


def concordance(Tg: pd.DataFrame, Pg: pd.DataFrame, t_members: dict, p_members: dict, cfg: dict) -> dict:
    """For each gene in both gene-level matrices: the correlation of its transcript with its protein across the samples."""
    method, fdr, min_n = cfg["method"], cfg["fdr"], cfg["min_samples"]
    common = sorted(set(Tg.columns) & set(Pg.columns))
    out = {"n_genes_t": int(Tg.shape[1]), "n_genes_p": int(Pg.shape[1]), "n_common": len(common), "rows": [],
           "n_tested": 0, "n_sig_pos": 0, "n_sig_neg": 0, "median_rho": None, "hist": None, "note": None}
    if not common:
        out["note"] = ("No gene is present in both the transcriptomics and the proteomics table (after identifier "
                       "resolution), so there is nothing to compare.")
        return out
    r, n = paired_correlation(Tg[common].to_numpy(), Pg[common].to_numpy(), method, min_n)
    p = pvalues_from_r(r, n)
    q = bh_qvalues(p, p)
    finite = np.isfinite(r)
    out["n_tested"] = int(finite.sum())
    order = np.argsort(np.where(np.isfinite(q), q, 2.0), kind="mergesort")
    order = sorted(order.tolist(), key=lambda i: (q[i] if np.isfinite(q[i]) else 2.0, -abs(r[i]) if np.isfinite(r[i]) else 0.0))
    for i in order:
        if not np.isfinite(r[i]):
            continue
        g = common[i]
        out["rows"].append({"Gene": g, "rho": float(r[i]), "p-value": float(p[i]), "FDR": float(q[i]), "n": int(n[i]),
                            "Transcript id(s)": ", ".join(t_members.get(g, [])), "Protein id(s)": ", ".join(p_members.get(g, []))})
    sig = np.isfinite(q) & (q <= fdr)
    out["n_sig_pos"] = int((sig & (r > 0)).sum())
    out["n_sig_neg"] = int((sig & (r < 0)).sum())
    if finite.any():
        out["median_rho"] = float(np.median(r[finite]))
        counts, edges = np.histogram(r[finite], bins=20, range=(-1.0, 1.0))
        out["hist"] = {"edges": [float(e) for e in edges], "counts": [int(c) for c in counts]}
    return out


# ---------------------------------------------------------------------------
# Cross-omics blocks
# ---------------------------------------------------------------------------

def cluster_order(M) -> list:
    """Leaf order of an average-linkage (UPGMA) clustering of the ROWS of M (Euclidean distance, NaN as 0). Own
    implementation: the app has no clustering code and scipy's is not used here."""
    X = np.nan_to_num(np.asarray(M, dtype=float))
    n = X.shape[0]
    if n <= 2:
        return list(range(n))
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2))
    clusters = {i: [i] for i in range(n)}
    nxt = n
    while len(clusters) > 1:
        ids = sorted(clusters)
        best = None
        for ai in range(len(ids)):
            for bi in range(ai + 1, len(ids)):
                a, b = ids[ai], ids[bi]
                d = float(np.mean([D[i, j] for i in clusters[a] for j in clusters[b]]))
                if best is None or d < best[0] - 1e-12:
                    best = (d, a, b)
        _, a, b = best
        clusters[nxt] = clusters.pop(a) + clusters.pop(b)
        nxt += 1
    return list(clusters[next(iter(clusters))])


def run_block(name_a: str, A: pd.DataFrame, name_b: str, B: pd.DataFrame, labels_a: dict, labels_b: dict, cfg: dict,
              progress=None) -> dict:
    """One omics pair: every feature of A against every feature of B."""
    method, fdr, rho_cut, min_n = cfg["method"], cfg["fdr"], cfg["rho_cut"], cfg["min_samples"]
    p_, q_ = A.shape[1], B.shape[1]
    if p_ * q_ > MAX_PAIRS_PER_BLOCK:
        raise ValueError(f"{OMICS_LABELS[name_a]} x {OMICS_LABELS[name_b]} would be {p_ * q_:,} pairs (limit "
                         f"{MAX_PAIRS_PER_BLOCK:,}): lower the number of top-variance features.")
    R, N = correlate_all(A.to_numpy(), B.to_numpy(), method, min_n, progress)
    P = pvalues_from_r(R, N)
    finite = np.isfinite(P)
    n_tests = int(finite.sum())
    pthr, n_sig = bh_threshold(P, fdr)
    mask = finite & (np.abs(R) >= rho_cut) & (P <= pthr) if n_sig else np.zeros(R.shape, dtype=bool)
    ii, jj = np.nonzero(mask)
    n_reported = int(ii.size)
    qv = bh_qvalues(P, P[ii, jj]) if n_reported else np.array([])
    order = sorted(range(n_reported), key=lambda k: (qv[k], -abs(float(R[ii[k], jj[k]]))))
    truncated = n_reported > MAX_REPORTED_ROWS
    order = order[:MAX_REPORTED_ROWS]
    ids_a, ids_b = list(A.columns), list(B.columns)
    rows = [{"Feature A": labels_a.get(ids_a[ii[k]], ids_a[ii[k]]), "Feature B": labels_b.get(ids_b[jj[k]], ids_b[jj[k]]),
             "Omics A": OMICS_LABELS[name_a], "Omics B": OMICS_LABELS[name_b], "rho": float(R[ii[k], jj[k]]),
             "p-value": float(P[ii[k], jj[k]]), "FDR": float(qv[k]), "n": int(N[ii[k], jj[k]]),
             "a_id": ids_a[ii[k]], "b_id": ids_b[jj[k]]} for k in order]
    out = {"omics_a": name_a, "omics_b": name_b, "n_a": p_, "n_b": q_, "n_tests": n_tests, "n_sig_fdr": int(n_sig),
           "n_reported": n_reported, "truncated": truncated, "expected_false": float(qv.sum()) if n_reported else 0.0,
           "n_too_few_samples": int(R.size - np.isfinite(R).sum()), "rows": rows, "heat": None, "edges": []}
    if n_reported:
        cnt_a, cnt_b = np.bincount(ii, minlength=p_), np.bincount(jj, minlength=q_)
        maxabs_a = np.zeros(p_)
        maxabs_b = np.zeros(q_)
        absr = np.abs(R[ii, jj])
        np.maximum.at(maxabs_a, ii, absr)
        np.maximum.at(maxabs_b, jj, absr)
        k_ = int(cfg["heat_k"])
        rows_sel = sorted(sorted(range(p_), key=lambda x: (-cnt_a[x], -maxabs_a[x], x))[:k_])
        cols_sel = sorted(sorted(range(q_), key=lambda x: (-cnt_b[x], -maxabs_b[x], x))[:k_])
        rows_sel = [x for x in rows_sel if cnt_a[x] > 0]
        cols_sel = [x for x in cols_sel if cnt_b[x] > 0]
        sub = R[np.ix_(rows_sel, cols_sel)].astype(float)
        ro, co = cluster_order(sub), cluster_order(sub.T)
        rows_sel, cols_sel = [rows_sel[i] for i in ro], [cols_sel[i] for i in co]
        sub = R[np.ix_(rows_sel, cols_sel)].astype(float)
        sig = mask[np.ix_(rows_sel, cols_sel)]
        out["heat"] = {"rows": [labels_a.get(ids_a[x], ids_a[x]) for x in rows_sel],
                       "cols": [labels_b.get(ids_b[x], ids_b[x]) for x in cols_sel],
                       "z": [[None if not np.isfinite(v) else round(float(v), 4) for v in r] for r in sub],
                       "sig": [[bool(v) for v in r] for r in sig]}
        top = sorted(range(n_reported), key=lambda k: -abs(float(R[ii[k], jj[k]])))[:int(cfg["net_k"])]
        out["edges"] = [{"a": ids_a[ii[k]], "b": ids_b[jj[k]], "a_label": labels_a.get(ids_a[ii[k]], ids_a[ii[k]]),
                         "b_label": labels_b.get(ids_b[jj[k]], ids_b[jj[k]]), "rho": float(R[ii[k], jj[k]]),
                         "fdr": float(qv[k])} for k in top]
    return out


def layered_layout(columns: "list[list]", edges: "list[tuple]", sweeps: int = 8) -> dict:
    """Place nodes in one column per omics. Rows are ordered by the barycentre of the neighbours (a few left-right
    sweeps) to cut edge crossings. columns: node ids per column, in column order; edges: (u, v) pairs.
    Returns {node: (column, row)}."""
    col_of = {n: c for c, nodes in enumerate(columns) for n in nodes}
    nbrs: dict = {n: [] for n in col_of}
    for u, v in edges:
        if u in col_of and v in col_of:
            nbrs[u].append(v)
            nbrs[v].append(u)
    order = [sorted(nodes, key=lambda n: (-len(nbrs[n]), str(n))) for nodes in columns]
    pos = {n: i for nodes in order for i, n in enumerate(nodes)}

    def reorder(c):
        size = max(len(order[c]) - 1, 1)
        def key(n):
            ps = [pos[m] / max(len(order[col_of[m]]) - 1, 1) for m in nbrs[n]]
            return (sum(ps) / len(ps) if ps else pos[n] / size, str(n))
        order[c] = sorted(order[c], key=key)
        for i, n in enumerate(order[c]):
            pos[n] = i
    for _ in range(sweeps):
        for c in list(range(len(columns))) + list(range(len(columns) - 2, -1, -1)):
            reorder(c)
    return {n: (col_of[n], pos[n]) for n in col_of}


def count_crossings(layout: dict, edges: "list[tuple]") -> int:
    """Pairs of edges that run between the same two columns and cross."""
    by_cols: dict = {}
    for u, v in edges:
        (cu, ru), (cv, rv) = layout[u], layout[v]
        if cu == cv:
            continue
        if cu > cv:
            (cu, ru), (cv, rv) = (cv, rv), (cu, ru)
        by_cols.setdefault((cu, cv), []).append((ru, rv))
    total = 0
    for es in by_cols.values():
        for i in range(len(es)):
            for j in range(i + 1, len(es)):
                if (es[i][0] - es[j][0]) * (es[i][1] - es[j][1]) < 0:
                    total += 1
    return total


def feature_labels(columns, gene_of: "dict | None") -> dict:
    """{feature id: label}: the gene symbol where the feature resolved to one, else the id itself; made unique."""
    labels, seen = {}, {}
    for c in columns:
        g = str((gene_of or {}).get(c, "") or "").strip()
        base = g if g else str(c)
        seen[base] = seen.get(base, 0) + 1
        labels[c] = base if seen[base] == 1 else f"{base} ({seen[base]})"
    return labels


def run_correlation(mats: dict, gene_of: dict, cfg: "dict | None" = None, progress=None) -> dict:
    """The whole analysis. mats: {omics: samples x features}, ALREADY aligned to the same samples (align_samples) and
    cleaned of unusable features; gene_of: {omics: {feature id: gene symbol}} for transcriptomics / proteomics.
    Returns {"config", "samples", "blocks", "concordance", "notes", "network"}; the matrices stay with the caller."""
    cfg = {**CONFIG_DEFAULTS, **(cfg or {})}
    _check_method(cfg["method"])
    omics = [o for o in OMICS_ORDER if o in mats]
    if len(omics) < 2:
        raise ValueError("Correlation needs at least two omics.")
    samples = list(mats[omics[0]].index)
    notes: list = []
    labels = {o: feature_labels(mats[o].columns, gene_of.get(o)) for o in omics}
    reduced, report = {}, {}
    for o in omics:
        reduced[o], report[o] = prepare_matrix(mats[o], False, MAX_MISSING_FRACTION, int(cfg["top_variance"]))
        if report[o]["n_variance_filtered"]:
            notes.append(f"{OMICS_LABELS[o]}: the {report[o]['n_out']} most variable of {report[o]['n_kept_before_variance']} "
                         "features are used for the cross-omics pairs (the same-gene test uses all genes).")
    blocks = {}
    pairs = [(a, b) for i, a in enumerate(omics) for b in omics[i + 1:]]
    for k, (a, b) in enumerate(pairs):
        # progress(fraction of the whole correlation step, text) - the blocks are the slow part
        step = None if progress is None else (lambda f, _k=k, _a=a, _b=b: progress(
            (_k + f) / len(pairs), f"{OMICS_LABELS[_a]} x {OMICS_LABELS[_b]}: {round(100 * f)}%"))
        blocks[f"{a}|{b}"] = run_block(a, reduced[a], b, reduced[b], labels[a], labels[b], cfg, step)
    conc = None
    if "transcriptomics" in omics and "proteomics" in omics:
        Tg, ti = collapse_by_gene(mats["transcriptomics"], gene_of.get("transcriptomics", {}))
        Pg, pi = collapse_by_gene(mats["proteomics"], gene_of.get("proteomics", {}))
        conc = concordance(Tg, Pg, ti["members"], pi["members"], cfg)
        conc["transcript_info"] = {k: v for k, v in ti.items() if k != "members"}
        conc["protein_info"] = {k: v for k, v in pi.items() if k != "members"}
        conc["_matrices"] = (Tg, Pg)
    # network: the top pairs of every block, one column per omics
    net_edges, nodes_by = [], {o: {} for o in omics}
    for key, blk in blocks.items():
        for e in blk["edges"]:
            ua, ub = f"{blk['omics_a']}:{e['a']}", f"{blk['omics_b']}:{e['b']}"
            nodes_by[blk["omics_a"]][ua] = e["a_label"]
            nodes_by[blk["omics_b"]][ub] = e["b_label"]
            net_edges.append({"u": ua, "v": ub, "rho": e["rho"], "fdr": e["fdr"], "block": key})
    cols = [[n for n in nodes_by[o]] for o in omics if nodes_by[o]]
    col_omics = [o for o in omics if nodes_by[o]]
    layout = layered_layout(cols, [(e["u"], e["v"]) for e in net_edges]) if net_edges else {}
    network = {"omics": col_omics, "edges": net_edges,
               "nodes": [{"id": n, "label": nodes_by[col_omics[c]][n], "omics": col_omics[c], "col": c, "row": r}
                         for n, (c, r) in layout.items()],
               "crossings": count_crossings(layout, [(e["u"], e["v"]) for e in net_edges]) if net_edges else 0}
    return {"config": cfg, "samples": samples, "omics": omics, "blocks": blocks, "concordance": conc, "notes": notes,
            "network": network, "feature_report": report, "labels": labels}


def pair_vectors(df_a: pd.DataFrame, fa: str, df_b: pd.DataFrame, fb: str) -> dict:
    """The two features' values per sample, for the scatter plot (None where missing)."""
    if fa not in df_a.columns or fb not in df_b.columns:
        raise KeyError("Unknown feature.")
    x, y = df_a[fa].to_numpy(dtype=float), df_b[fb].to_numpy(dtype=float)
    f = lambda v: None if not np.isfinite(v) else float(v)
    return {"samples": [str(s) for s in df_a.index], "x": [f(v) for v in x], "y": [f(v) for v in y]}
