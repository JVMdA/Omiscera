# -*- coding: utf-8 -*-
"""Unit tests for Multiomics_correlation (Correlation-based Integration engine), offline.

THE POINT OF THIS FILE is that every number is checked against something INDEPENDENT of the engine:
  - ranks, correlations (with ties and missing values) and p-values against a pure-python reference written here,
    whose p-value comes from a hand-written regularised incomplete beta (a continued fraction), not from scipy;
  - Benjamini-Hochberg against Genomics_module's own list-based implementation;
  - the clustering order against a naive average-linkage written here;
  - and, most importantly, NULL CALIBRATION: on matrices of independent noise the p-values must be uniform and the
    share below 0.05 must be about 5%, and BH must find (almost) nothing. That is the property that makes a
    correlation method honest rather than merely plausible.
Rules (each mutation-tested): samples matched by id and only the shared ones used, reported when dropped; features
with too many missing values or constant dropped and counted; pairwise-complete samples with each pair's own n;
BH per block over EVERY tested pair; a pair is reported only if q <= FDR AND |rho| >= cut; same-gene concordance
uses all genes while the blocks use the top-variance features; duplicates collapse by mean; no matrix product anywhere.
"""
import math
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm  # noqa: E402
import Multiomics_correlation as mc  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def raises(fn, *a, **k):
    try:
        fn(*a, **k)
    except Exception as exc:  # noqa: BLE001
        return type(exc).__name__, str(exc)
    return None, ""


# ------------------------------------------------------------------ independent references
def ref_rank(xs):
    idx = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and xs[idx[j + 1]] == xs[idx[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[idx[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def ref_pearson(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    n = len(pairs)
    if n < 2:
        return float("nan"), n
    mx, my = sum(a for a, _ in pairs) / n, sum(b for _, b in pairs) / n
    sxx = sum((a - mx) ** 2 for a, _ in pairs)
    syy = sum((b - my) ** 2 for _, b in pairs)
    sxy = sum((a - mx) * (b - my) for a, b in pairs)
    if sxx <= 0 or syy <= 0:
        return float("nan"), n
    return max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy))), n


def ref_spearman(x, y):
    """The standard pairwise-complete Spearman (R's cor(use = "pairwise.complete.obs", method = "spearman"), pandas):
    keep the samples BOTH have, rank each vector over those samples, Pearson of the ranks."""
    pairs = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    if not pairs:
        return float("nan"), 0
    return ref_pearson(ref_rank([a for a, _ in pairs]), ref_rank([b for _, b in pairs]))


def old_shortcut_spearman(x, y):
    """What the engine USED to do (each vector ranked over its own values, then Pearson over the shared samples) - kept
    only to show the tests can tell the two apart."""
    def rk(v):
        ok = [i for i, a in enumerate(v) if a == a]
        out = [float("nan")] * len(v)
        for i, r_ in zip(ok, ref_rank([v[i] for i in ok])):
            out[i] = r_
        return out
    return ref_pearson(rk(x), rk(y))


def _betacf(a, b, x):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h


def ref_betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def ref_p(r, n):
    if r != r or n < 3:
        return float("nan")
    if abs(r) >= 1.0:
        return 0.0
    df = n - 2
    t = r * math.sqrt(df / (1.0 - r * r))
    return ref_betai(df / 2.0, 0.5, df / (df + t * t))


def close(a, b, rel=1e-7, ab=1e-12):
    if a != a and b != b:
        return True
    return abs(a - b) <= max(ab, rel * max(abs(a), abs(b)))


CSV_ROWS_SAMPLES = (b"Sample,Group,Citrate,Lactate,Pyruvate\n"
                    b"Control_1,Control,1.0,2.0,3.0\nControl_2,Control,1.5,2.5,3.2\nControl_3,Control,0.9,2.2,3.9\n"
                    b"Treated_1,Treated,2.0,1.0,2.0\nTreated_2,Treated,2.2,1.1,2.4\n")
CSV_COLS_SAMPLES = (b"Gene,S1,S2,S3,S4,S5,S6\n" + b"\n".join(
    (f"G{i}," + ",".join(str(round(math.sin(i * 7 + j * 3) + i, 4)) for j in range(6))).encode() for i in range(10)) + b"\n")

print("=== reading a matrix ===")
r1 = mc.parse_matrix(CSV_ROWS_SAMPLES)
df1, m1 = r1["df"], r1["meta"]
check(df1.shape == (5, 3) and list(df1.columns) == ["Citrate", "Lactate", "Pyruvate"] and list(df1.index)[0] == "Control_1"
      and m1["orientation"] == "rows" and m1["set_aside_columns"] == ["Group"],
      f"samples in rows with a Group column: the Group column is set aside, not data; auto says rows ({m1['orientation_reason']})")
check(m1["group_columns"] == {"Group": {"Control_1": "Control", "Control_2": "Control", "Control_3": "Control", "Treated_1": "Treated", "Treated_2": "Treated"}},
      "the set-aside Group column is kept as per-sample labels (MOFA colours and tests its factors by them)")
check(mc.parse_matrix(CSV_COLS_SAMPLES)["meta"]["group_columns"] == {}, "an expression-style table (samples in columns) has no per-sample label columns")
_lines = CSV_ROWS_SAMPLES.decode().splitlines()
_uniq = "\n".join([_lines[0]] + [l.replace(",Control,", f",G{i},").replace(",Treated,", f",G{i},") for i, l in enumerate(_lines[1:])]) + "\n"
check(mc.parse_matrix(_uniq.encode())["meta"]["group_columns"] == {}, "SENSITIVITY: a label column whose every value is different is not offered as a group")
_one = "\n".join([_lines[0]] + [l.replace(",Treated,", ",Control,") for l in _lines[1:]]) + "\n"
check(mc.parse_matrix(_one.encode())["meta"]["group_columns"] == {}, "SENSITIVITY: a label column with a single value is not offered as a group")
_gl = CSV_COLS_SAMPLES.decode().splitlines()
_withtype = "\n".join(["Gene,Type," + _gl[0].split(",", 1)[1]] + [f"{l.split(chr(44), 1)[0]},{chr(97) if i % 2 else chr(98)}," + l.split(",", 1)[1] for i, l in enumerate(_gl[1:])]) + "\n"
_mt = mc.parse_matrix(_withtype.encode())["meta"]
check(_mt["orientation"] == "columns" and "Type" in _mt["set_aside_columns"] and _mt["group_columns"] == {},
      "SENSITIVITY: a text column in a samples-in-columns table describes FEATURES, so it is not a sample grouping")
r2 = mc.parse_matrix(CSV_COLS_SAMPLES)
check(r2["df"].shape == (6, 10) and r2["meta"]["orientation"] == "columns" and list(r2["df"].index) == [f"S{j}" for j in range(1, 7)],
      "an expression-style table (many rows, few columns) is read as samples in the columns and transposed to samples x features")
r3 = mc.parse_matrix(CSV_COLS_SAMPLES, orientation="rows")
check(r3["df"].shape == (10, 6) and list(r3["df"].index)[0] == "G0", "SENSITIVITY: the orientation you set overrides the guess")
check(mc.parse_matrix(CSV_ROWS_SAMPLES, "columns")["df"].shape == (3, 5), "...and an explicit 'columns' transposes")
check(sorted(r1["meta"]["alt_samples"]) == sorted(["Citrate", "Lactate", "Pyruvate"]), "the labels of the OTHER axis are kept, so a wrong orientation can be diagnosed")
dup = b"Sample,A,B\nS1,1,2\nS1,3,4\nS2,5,6\n"
name, msg = raises(mc.parse_matrix, dup, "rows")
check(name == "ValueError" and "Duplicate sample ids" in msg, f"duplicate sample ids are an error that says what to try ({msg[:80]})")
dupf = mc.parse_matrix(b"Gene,S1,S2,S3\nA,1,2,3\nA,2,3,5\nB,3,1,2\nC,5,5,1\nD,1,4,2\n")
check(list(dupf["df"].columns) == ["A", "A (2)", "B", "C", "D"] and dupf["meta"]["n_repeated_feature_ids"] == 1,
      "a repeated FEATURE id is kept as a separate feature, renamed, and counted")
check(dupf["meta"]["original_ids"] == {"A (2)": "A"}, "the renamed feature remembers its original id (so it can still be resolved to a gene)")


def _orient(header_extra, values_extra, n_rows=40, n_cols=6, extra_first=False):
    rows = []
    for i in range(n_rows):
        nums = ",".join(str(round(math.sin(i * 3 + j) + 5, 3)) for j in range(n_cols))
        rows.append(f"F{i},{values_extra(i)},{nums}" if extra_first else f"F{i},{nums},{values_extra(i)}")
    head = ",".join(f"S{j}" for j in range(n_cols))
    header = f"Id,{header_extra},{head}" if extra_first else f"Id,{head},{header_extra}"
    return mc.parse_matrix((header + "\n" + "\n".join(rows) + "\n").encode())["meta"]


for col, vals, want, why in [
        ("gene_biotype", lambda i: ["protein_coding", "lncRNA"][i % 2], "columns", "an RNA table's gene_biotype annotation"),
        ("Class", lambda i: ["Lipid", "Amino acid", "Sugar"][i % 3], "columns", "a compound table's chemical Class"),
        ("Description", lambda i: f"protein {i % 4}", "columns", "a free-text Description"),
        ("Protein names", lambda i: ["Kinase", "Ligase"][i % 2], "columns", "a 'Protein names' annotation"),
        ("Compound group", lambda i: ["Lipids", "Acids"][i % 2], "columns", "a 'Compound group' annotation (a group word ON a feature)")]:
    m_ = _orient(col, vals)
    check(m_["orientation"] == "columns" and col in m_["set_aside_columns"],
          f"SENSITIVITY: {why} is set aside and does NOT turn 40 features into samples (got {m_['orientation']})")
for col, vals, why in [("Label", lambda i: ["Control", "Case"][i % 2], "MetaboAnalyst's 'Label' column"),
                       ("Sample type", lambda i: ["QC", "Study"][i % 2], "a 'Sample type' column"),
                       ("Treatment", lambda i: ["A", "B", "C"][i % 3], "a Treatment column")]:
    m_ = _orient(col, vals, extra_first=True)
    check(m_["orientation"] == "rows" and col in m_["orientation_reason"],
          f"{why} labels each row as a sample, even with more rows than columns (reason names it: {m_['orientation_reason'][:60]})")
m_ = _orient("Group", lambda i: f"g{i}", extra_first=True)
check(m_["orientation"] == "columns", "a 'Group' column whose every value is different is not a sample grouping (nothing repeats)")
check(raises(mc.parse_matrix, b"")[0] == "ValueError" and raises(mc.parse_matrix, b"a,b\n1,2\n")[0] == "ValueError"
      and raises(mc.parse_matrix, b"id,x,y\na,hello,world\nb,foo,bar\nc,x,y\n")[0] == "ValueError", "an empty, a too-small or a non-numeric file is refused")
tab = mc.parse_matrix(b"Gene\tS1\tS2\tS3\nA\t1\t2\t3\nB\t2\t3\t5\nC\t3\t1\t2\nD\t5\t5\t1\n")
check(tab["df"].shape == (3, 4), "a tab-separated file is read too")
blank = mc.parse_matrix(b"Gene,S1,S2,S3\nA,1,2,3\nB,,,\nC,3,1,2\nD,5,5,1\n")
check(list(blank["df"].columns) == ["A", "C", "D"], "an all-empty feature row is dropped")
mixed = mc.parse_matrix(b"Gene,S1,S2,S3\nA,1,2,3\nB,2,x,5\nC,3,1,2\nD,5,5,1\n")
check(np.isnan(mixed["df"].loc["S2", "B"]) and mixed["df"].loc["S1", "B"] == 2.0, "a non-numeric cell inside a numeric column is a missing value")

print("=== cleaning ===")
rng = np.random.RandomState(0)
raw = pd.DataFrame(rng.rand(10, 6), index=[f"S{i}" for i in range(10)], columns=list("ABCDEF"))
raw["C"] = 3.0                                     # constant
raw.loc[raw.index[:4], "D"] = np.nan               # 40% missing
raw.loc[raw.index[:1], "E"] = np.nan               # 10% missing: kept
clean, rep = mc.prepare_matrix(raw)
check(list(clean.columns) == ["A", "B", "E", "F"] and rep["n_constant"] == 1 and rep["n_too_many_missing"] == 1,
      f"SENSITIVITY: a constant feature and one with 40% missing are dropped and counted; 10% missing is kept ({list(clean.columns)})")
c2, r2_ = mc.prepare_matrix(raw, top_variance=2)
v = raw[["A", "B", "E", "F"]].var()
check(list(c2.columns) == sorted(v.sort_values(ascending=False).index[:2], key=list("ABEF").index) and r2_["n_variance_filtered"] == 2,
      f"the top-variance features are kept, in their original order ({list(c2.columns)})")
c3, _ = mc.prepare_matrix(raw[["A", "B"]] * 100, log2=True)
check(np.allclose(c3["A"], np.log2(raw["A"] * 100 + 1)), "log2(x + 1) is applied on request")
check(raises(mc.prepare_matrix, raw[["A", "B"]] - 0.5, True)[0] == "ValueError", "log2 of negative values is refused (the data are already log-scaled)")

print("=== matching samples ===")
a = pd.DataFrame(np.arange(12.0).reshape(4, 3), index=["S1", "S2 ", "s3", "S4"], columns=list("xyz"))
b = pd.DataFrame(np.arange(8.0).reshape(4, 2), index=["S4", "S2", "S1", "S9"], columns=list("uv"))
al = mc.align_samples({"transcriptomics": a, "proteomics": b}, min_samples=2)
check(al["samples"] == ["S1", "S2", "S4"] and list(al["mats"]["proteomics"].index) == ["S1", "S2", "S4"]
      and np.allclose(al["mats"]["proteomics"].loc["S1"], b.loc["S1"]) and np.allclose(al["mats"]["transcriptomics"].loc["S2"], a.loc["S2 "]),
      "samples are matched by id with whitespace ignored, in the first omics' order, and each row keeps ITS OWN values")
check(al["report"]["transcriptomics"]["dropped"] == ["s3"] and al["report"]["proteomics"]["dropped"] == ["S9"],
      "SENSITIVITY: the dropped samples are named per omics (never silently)")
alc = mc.align_samples({"transcriptomics": a, "proteomics": b}, case_insensitive=True, min_samples=2)
check(alc["samples"] == ["S1", "S2", "S4"], "case-insensitive matching changes nothing when no ids differ only by case")
b2 = b.copy(); b2.index = ["S4", "S2", "S3", "S9"]       # the RNA table calls it 's3', this one 'S3'
alc_off = mc.align_samples({"transcriptomics": a, "proteomics": b2}, case_insensitive=False, min_samples=2)
alc_on = mc.align_samples({"transcriptomics": a, "proteomics": b2}, case_insensitive=True, min_samples=2)
check(len(alc_off["samples"]) == 2 and "s3" not in [s.lower() for s in alc_off["samples"]] and alc_off["report"]["transcriptomics"]["dropped"] == ["S1", "s3"],
      "SENSITIVITY: with case respected, 's3' and 'S3' are different samples (and the ones left out are named)")
check(len(alc_on["samples"]) == 3 and "s3" in [s.lower() for s in alc_on["samples"]] and alc_on["report"]["transcriptomics"]["dropped"] == ["S1"],
      "SENSITIVITY: with case ignored, 's3' matches 'S3'")
name, msg = raises(mc.align_samples, {"transcriptomics": a, "proteomics": b}, False, 5)
check(name == "ValueError" and "at least 5" in msg, f"fewer shared samples than the minimum is refused ({msg[:70]})")
tr = pd.DataFrame(np.arange(18.0).reshape(6, 3), index=[f"S{i}" for i in range(6)], columns=list("abc"))
wrong = pd.DataFrame(np.arange(18.0).reshape(3, 6), index=list("abc"), columns=[f"S{i}" for i in range(6)])
wrong.attrs["alt_samples"] = [f"S{i}" for i in range(6)]
name, msg = raises(mc.align_samples, {"transcriptomics": tr, "proteomics": wrong}, False, 5)
check(name == "ValueError" and "other orientation" in msg, "SENSITIVITY: when the other orientation of a file would match, the error says so")
check(len(mc.align_samples({"transcriptomics": tr, "proteomics": tr.copy()}, min_samples=5)["warnings"]) == 1, "fewer than 10 shared samples warns that the p-values are approximate")
dupk = a.copy(); dupk.index = ["S1", "S1 ", "S3", "S4"]
check(raises(mc.align_samples, {"transcriptomics": dupk, "proteomics": b}, False, 2)[0] == "ValueError", "two samples with the same id after trimming are refused")

print("=== correlation equals an independent reference ===")
rs = np.random.RandomState(5)
n = 14
A = rs.randn(n, 9)
B = rs.randn(n, 7)
A[:, 1] = np.round(A[:, 1])                       # many ties
B[:, 2] = np.round(B[:, 2] * 2)
A[:, 3] = B[:, 0] * 2 + 1                          # exact linear: r = 1
A[:, 4] = -B[:, 1]                                 # exact: r = -1
for method, ref in (("pearson", ref_pearson), ("spearman", ref_spearman)):
    R, N = mc.correlate_all(A, B, method, 5)
    worst = 0.0
    for i in range(A.shape[1]):
        for j in range(B.shape[1]):
            rr, nn = ref(list(A[:, i]), list(B[:, j]))
            worst = max(worst, abs(rr - float(R[i, j])) if rr == rr else (0 if np.isnan(R[i, j]) else 1))
    check(worst < 2e-6 and (N == n).all(), f"{method}: all {A.shape[1] * B.shape[1]} correlations equal an independent computation (worst difference {worst:.1e}), no missing values")
R, N = mc.correlate_all(A, B, "pearson", 5)
check(abs(R[3, 0] - 1) < 1e-6 and abs(R[4, 1] + 1) < 1e-6, "exact linear relations give r = +1 and -1")
A_m, B_m = A.copy(), B.copy()
A_m[[0, 3], 2] = np.nan
B_m[[3, 5, 6], 0] = np.nan
B_m[:, 5] = np.nan
B_m[:10, 6] = np.nan
for method, ref in (("pearson", ref_pearson), ("spearman", ref_spearman)):
    R, N = mc.correlate_all(A_m, B_m, method, 5)
    worst, nbad = 0.0, 0
    for i in range(A_m.shape[1]):
        for j in range(B_m.shape[1]):
            rr, nn = ref(list(A_m[:, i]), list(B_m[:, j]))
            if nn < 5:
                nbad += 0 if np.isnan(R[i, j]) else 1
                continue
            # NaN-safe: max(x, nan) silently returns x, so a value that became NaN must count as wrong explicitly
            worst = max(worst, abs(rr - float(R[i, j])) if np.isfinite(R[i, j]) else 1.0)
            nbad += 0 if int(N[i, j]) == nn else 1
    check(worst < 2e-6 and nbad == 0,
          f"SENSITIVITY: {method} with missing values - each pair uses ITS OWN complete samples (n differs by pair) and equals the reference (worst {worst:.1e}, {nbad} wrong n)")
gap = max(abs(old_shortcut_spearman(list(A_m[:, i]), list(B_m[:, j]))[0] - ref_spearman(list(A_m[:, i]), list(B_m[:, j]))[0])
          for i in range(A_m.shape[1]) for j in range(5))
check(gap > 0.01, f"SENSITIVITY: the check above can tell exact pairwise Spearman from ranking each feature on its own (they differ by up to {gap:.3f} here)")
# missing values AND ties AND both sides missing in different samples, at a realistic size
rt = np.random.RandomState(11)
At = np.round(rt.randn(12, 25) * 2)                # heavy ties
Bt = np.round(rt.randn(12, 30) * 2)
At[rt.rand(12, 25) < 0.12] = np.nan
Bt[rt.rand(12, 30) < 0.12] = np.nan
Rt, Nt = mc.correlate_all(At, Bt, "spearman", 5)
worst_t = 0.0
for i in range(25):
    for j in range(30):
        rr_, nn_ = ref_spearman(list(At[:, i]), list(Bt[:, j]))
        if nn_ >= 5 and rr_ == rr_:
            worst_t = max(worst_t, abs(rr_ - float(Rt[i, j])) if np.isfinite(Rt[i, j]) else 1.0)
            worst_t = max(worst_t, 0 if int(Nt[i, j]) == nn_ else 1)
check(worst_t < 2e-9, f"SENSITIVITY: tie-heavy data with ~12% missing on BOTH sides - all 750 Spearman values equal the per-pair reference (worst {worst_t:.1e})")
rp, np_ = mc.paired_correlation(At[:, :25], Bt[:, :25], "spearman", 5)
check(all((abs(rp[k] - ref_spearman(list(At[:, k]), list(Bt[:, k]))[0]) < 2e-9) or np.isnan(rp[k]) for k in range(25)),
      "SENSITIVITY: the same-gene (paired) Spearman is ranked over each pair's shared samples too")
Xr = At.copy()
rk_ = mc._rank_columns(Xr)
check(all(np.allclose(rk_[~np.isnan(Xr[:, j]), j], ref_rank(list(Xr[~np.isnan(Xr[:, j]), j]))) and np.isnan(rk_[np.isnan(Xr[:, j]), j]).all()
          for j in range(Xr.shape[1])), "the vectorised column ranking equals the reference (ties averaged, missing stays missing)")
check(int(N[2, 0]) == n - 4 and int(N[2, 1]) == n - 2 and np.isnan(R[:, 5]).all() and np.isnan(R[:, 6]).all() and (N[:, 6] <= 4).all(),
      "the pair counts reflect the missing cells; a column with too few common samples is NaN, never a number from 4 points")
check(np.isnan(mc.correlate_all(np.c_[A[:, 0], np.ones(n)], B[:, :2], "pearson", 5)[0][1]).all(), "a constant feature gives NaN, not 0 or an error")
r_pair, n_pair = mc.paired_correlation(A[:, :7], B, "spearman", 5)
check(all(abs(r_pair[k] - ref_spearman(list(A[:, k]), list(B[:, k]))[0]) < 2e-6 for k in range(7)) and (n_pair == n).all(),
      "paired (same-column) correlation equals the reference")
check(raises(mc.correlate_all, A, B, "kendall")[0] == "ValueError", "an unknown method is refused")
check(all(abs(a_ - b_) < 1e-12 for a_, b_ in zip(mc._rankdata(np.array([3.0, 1.0, 1.0, 2.0, 3.0])), ref_rank([3.0, 1.0, 1.0, 2.0, 3.0]))),
      "average ranks equal the reference with ties")

print("=== p-values and Benjamini-Hochberg ===")
rr = np.array([0.7, -0.7, 0.0, 0.3, 0.95, 0.999])
nn = np.array([12, 12, 12, 6, 30, 5])
pp = mc.pvalues_from_r(rr, nn)
check(all(close(float(pp[k]), ref_p(float(rr[k]), int(nn[k])), 1e-9) for k in range(6)),
      "SENSITIVITY: the p-values equal an independent computation from a hand-written incomplete beta (r = 0.7, n = 12 -> 0.01126)")
check(abs(float(pp[0]) - 0.011257326) < 1e-8 and pp[0] == pp[1] and abs(float(pp[2]) - 1) < 1e-12, "r = 0 gives p = 1, and sign does not matter")
pe = mc.pvalues_from_r(np.array([1.0, -1.0, np.nan, 0.5]), np.array([10, 10, 10, 2]))
check(pe[0] == 0.0 and pe[1] == 0.0 and np.isnan(pe[2]) and np.isnan(pe[3]), "r = +-1 gives p = 0, a missing r or n < 3 gives NaN")
pr = np.random.RandomState(2).rand(500)
pr[::7] = np.round(pr[::7], 2)                         # ties
q = mc.bh_qvalues(pr, pr)
ref_q = np.array(gm.benjamini_hochberg(list(pr)))
check(np.allclose(q, ref_q, atol=1e-12), "BH q-values equal Genomics_module's own implementation (ties included)")
pf = np.r_[pr, np.nan, np.nan]
check(np.allclose(mc.bh_qvalues(pf, pr), ref_q), "NaN p-values are not part of the family (m counts the tested pairs only)")
check(np.allclose(mc.bh_qvalues(pr, pr[:10]), ref_q[:10]), "q-values of a SUBSET are computed inside the full family, not the subset")
for fdr in (0.01, 0.05, 0.2):
    thr, k = mc.bh_threshold(pr, fdr)
    check(k == int((ref_q <= fdr).sum()) and (k == 0 or abs(thr - float(np.max(pr[ref_q <= fdr]))) < 1e-15), f"the BH threshold at FDR {fdr} selects exactly the pairs with q <= {fdr} ({k})")
check(mc.bh_threshold(np.array([0.9, 0.8, 0.7]), 0.05) == (0.0, 0), "nothing significant -> (0, 0)")

print("=== NULL CALIBRATION: independent noise must look like noise ===")
for method in ("spearman", "pearson"):
    rn = np.random.RandomState(11)
    A0, B0 = rn.randn(20, 300), rn.randn(20, 300)
    R0, N0 = mc.correlate_all(A0, B0, method, 5)
    P0 = mc.pvalues_from_r(R0, N0).ravel()
    frac = float((P0 < 0.05).mean())
    ks = float(np.max(np.abs(np.sort(P0) - (np.arange(1, P0.size + 1) / P0.size))))
    thr, k = mc.bh_threshold(P0, 0.05)
    check(0.04 < frac < 0.06 and ks < 0.02,
          f"{method}, n = 20, 90,000 independent pairs: {frac:.1%} have p < 0.05 and the p-values are uniform (KS distance {ks:.3f}); BH finds {k}")
    check(k == 0, f"{method}: BH at FDR 5% reports NOTHING on pure noise")
rn = np.random.RandomState(12)
A0, B0 = rn.randn(12, 300), rn.randn(12, 300)
R0, N0 = mc.correlate_all(A0, B0, "spearman", 5)
P0 = mc.pvalues_from_r(R0, N0).ravel()
check(0.035 < float((P0 < 0.05).mean()) < 0.065, f"spearman at n = 12 is still calibrated within sampling noise ({float((P0 < 0.05).mean()):.1%} below 0.05)")

print("=== blocks and planted structure ===")
rp = np.random.RandomState(21)
n = 16
base = rp.randn(n)
idx = [f"S{i}" for i in range(n)]
T = pd.DataFrame(rp.randn(n, 300), index=idx, columns=[f"t{i}" for i in range(300)])
Pm = pd.DataFrame(rp.randn(n, 250), index=idx, columns=[f"p{i}" for i in range(250)])
Mm = pd.DataFrame(rp.randn(n, 40), index=idx, columns=[f"m{i}" for i in range(40)])
for k in range(3):
    T[f"t{k}"] = base * (1 + k) + 0.3 * rp.randn(n)
    Pm[f"p{k}"] = base * (1 + k) + 0.4 * rp.randn(n)
Mm["m0"] = base + 0.3 * rp.randn(n)
gene_of = {"transcriptomics": {c: f"GENE{c[1:]}" for c in T.columns}, "proteomics": {c: f"GENE{c[1:]}" for c in Pm.columns}}
res = mc.run_correlation({"transcriptomics": T, "proteomics": Pm, "metabolomics": Mm}, gene_of, {"top_variance": 1000})
b_tp = res["blocks"]["transcriptomics|proteomics"]
check(set(res["blocks"]) == {"transcriptomics|proteomics", "transcriptomics|metabolomics", "proteomics|metabolomics"}, "three omics give three blocks (one family each)")
check(b_tp["n_tests"] == 300 * 250 and res["blocks"]["transcriptomics|metabolomics"]["n_tests"] == 300 * 40,
      "SENSITIVITY: the BH family of a block is EVERY tested pair (300 x 250), not the pairs that were reported")
top = {(r_["Feature A"], r_["Feature B"]) for r_ in b_tp["rows"]}
check({("GENE0", "GENE0"), ("GENE1", "GENE1"), ("GENE2", "GENE2")} <= top, "the planted transcript-protein pairs are all reported")
check(all(abs(r_["rho"]) >= 0.6 and r_["FDR"] <= 0.05 for r_ in b_tp["rows"]), "every reported pair has |rho| >= the cut AND q <= the FDR")
check(any(r_["Feature A"] == "GENE0" and r_["Feature B"] == "m0" for r_ in res["blocks"]["transcriptomics|metabolomics"]["rows"]),
      "a planted gene-metabolite pair is found in its own block (labels use the gene symbol, metabolites keep their name)")
qs = [r_["FDR"] for r_ in b_tp["rows"]]
check(qs == sorted(qs) or all(qs[i] <= qs[i + 1] + 1e-12 for i in range(len(qs) - 1)), "pairs are listed by FDR")
ref_all = []
for i in range(300):
    for j in range(250):
        rr_, nn_ = ref_spearman(list(T.iloc[:, i]), list(Pm.iloc[:, j]))
        ref_all.append(ref_p(rr_, nn_))
ref_q_all = np.array(gm.benjamini_hochberg(ref_all)).reshape(300, 250)
reported = {(r_["a_id"], r_["b_id"]): r_["FDR"] for r_ in b_tp["rows"]}
ok = all(abs(reported[(f"t{i}", f"p{j}")] - ref_q_all[i, j]) < 1e-9 for (i, j) in [(0, 0), (1, 1), (2, 2)] if (f"t{i}", f"p{j}") in reported)
check(ok, "INDEPENDENT: the reported q-values equal BH over all 75,000 pairs computed by the pure-python reference")
check(abs(b_tp["expected_false"] - sum(qs)) < 1e-9 and b_tp["n_reported"] == len(b_tp["rows"]) and not b_tp["truncated"],
      "the expected number of false pairs is the sum of the reported q-values")
check(b_tp["heat"] is not None and len(b_tp["heat"]["rows"]) <= 40 and len(b_tp["heat"]["z"]) == len(b_tp["heat"]["rows"]) and len(b_tp["heat"]["z"][0]) == len(b_tp["heat"]["cols"]),
      "a heatmap block is returned, sized to its rows and columns")
check({r_ for r_ in b_tp["heat"]["rows"]} >= {"GENE0", "GENE1", "GENE2"}, "the heatmap contains the features with reported partners")
check(any(any(row) for row in b_tp["heat"]["sig"]) and all(isinstance(v, bool) for v in b_tp["heat"]["sig"][0]), "...and marks which cells are reported")
cut_hi = mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"rho_cut": 0.99, "top_variance": 1000})
check(len(cut_hi["blocks"]["transcriptomics|proteomics"]["rows"]) < len(b_tp["rows"]), "SENSITIVITY: a higher |rho| cutoff reports fewer pairs")
fdr_lo = mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"fdr": 1e-12, "rho_cut": 0.0, "top_variance": 1000})
check(len(fdr_lo["blocks"]["transcriptomics|proteomics"]["rows"]) < len(mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"fdr": 0.2, "rho_cut": 0.0, "top_variance": 1000})["blocks"]["transcriptomics|proteomics"]["rows"]),
      "SENSITIVITY: a stricter FDR reports fewer pairs")
nulls = mc.run_correlation({"transcriptomics": pd.DataFrame(np.random.RandomState(3).randn(16, 200), index=idx, columns=[f"t{i}" for i in range(200)]),
                            "proteomics": pd.DataFrame(np.random.RandomState(4).randn(16, 200), index=idx, columns=[f"p{i}" for i in range(200)])},
                           {}, {"rho_cut": 0.0})
check(nulls["blocks"]["transcriptomics|proteomics"]["n_reported"] == 0 and nulls["blocks"]["transcriptomics|proteomics"]["rows"] == [],
      "the whole pipeline on independent noise reports no pair")
tv = mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"top_variance": 50})
check(tv["blocks"]["transcriptomics|proteomics"]["n_a"] == 50 and tv["blocks"]["transcriptomics|proteomics"]["n_tests"] == 50 * 50 and tv["notes"],
      "the top-variance filter limits the cross-omics block and says so")
old = mc.MAX_PAIRS_PER_BLOCK
mc.MAX_PAIRS_PER_BLOCK = 1000
check(raises(mc.run_correlation, {"transcriptomics": T, "proteomics": Pm}, gene_of, {"top_variance": 1000})[0] == "ValueError", "a block of too many pairs is refused with advice")
mc.MAX_PAIRS_PER_BLOCK = old
old_r = mc.MAX_REPORTED_ROWS
mc.MAX_REPORTED_ROWS = 2
cap = mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"top_variance": 1000})["blocks"]["transcriptomics|proteomics"]
check(cap["truncated"] and len(cap["rows"]) == 2 and cap["n_reported"] > 2, "SENSITIVITY: more pairs than the row cap are counted but only the cap is listed, and it says it was cut")
mc.MAX_REPORTED_ROWS = old_r
Tm = T.copy()
Tm.iloc[:3, 3] = np.nan
resm = mc.run_correlation({"transcriptomics": Tm, "proteomics": Pm}, gene_of, {"top_variance": 1000})
check(resm["blocks"]["transcriptomics|proteomics"]["n_tests"] == 300 * 250, "missing values below the limit keep the feature in the family")
check(raises(mc.run_correlation, {"transcriptomics": T}, gene_of)[0] == "ValueError", "one omics is not a correlation")

print("=== same-gene concordance ===")
rc = res["concordance"]
check(rc["n_common"] == 250 and rc["n_tested"] == 250 and rc["n_genes_t"] == 300 and rc["n_genes_p"] == 250, "genes present in both tables are compared (250 of 300 transcripts have a protein)")
row0 = {r_["Gene"]: r_ for r_ in rc["rows"]}
check(all(row0[f"GENE{k}"]["rho"] > 0.8 and row0[f"GENE{k}"]["FDR"] < 0.01 for k in range(3)), "the planted concordant genes have high rho and small FDR")
g_all = [ref_spearman(list(T[f"t{k}"]), list(Pm[f"p{k}"]))[0] for k in range(3)]
check(all(abs(row0[f"GENE{k}"]["rho"] - g_all[k]) < 2e-6 for k in range(3)), "INDEPENDENT: their rho equals the reference")
check(sum(rc["hist"]["counts"]) == rc["n_tested"] and len(rc["hist"]["edges"]) == 21, "the histogram counts every tested gene")
check(abs(rc["median_rho"] - float(np.median([r_["rho"] for r_ in rc["rows"]]))) < 1e-12, "the median rho is the median of the gene correlations")
ref_pg = [ref_p(row0[g]["rho"], row0[g]["n"]) for g in row0]
q_ref = dict(zip(row0, gm.benjamini_hochberg(ref_pg)))
check(all(abs(row0[g]["FDR"] - q_ref[g]) < 1e-9 for g in row0), "INDEPENDENT: the concordance FDR is BH over the 250 genes only (its own family)")
check(rc["n_sig_pos"] >= 3 and rc["n_sig_neg"] == 0, f"significant positive / negative counts ({rc['n_sig_pos']} / {rc['n_sig_neg']})")
dg = pd.DataFrame({"a": [1.0, 2, 3, 4, 5, 6], "b": [2.0, 3, 5, 7, 8, 9], "c": [6.0, 5, 4, 3, 2, 1], "d": [1.0, 3, 2, 5, 4, 6]}, index=[f"S{i}" for i in range(6)])
coll, info = mc.collapse_by_gene(dg, {"a": "geneX", "b": "GENEX", "c": "", "d": "geneY"})
check(list(coll.columns) == ["GENEX", "GENEY"] and np.allclose(coll["GENEX"], (dg["a"] + dg["b"]) / 2) and info["n_unmapped"] == 1 and info["n_collapsed_genes"] == 1,
      "SENSITIVITY: features that resolve to one gene (any case) are averaged; a feature with no gene is left out and counted")
check(info["members"]["GENEX"] == ["a", "b"], "the members of a collapsed gene are listed")
dgn = dg.copy(); dgn.loc["S0", "a"] = np.nan
check(abs(mc.collapse_by_gene(dgn, {"a": "G", "b": "G"})[0].loc["S0", "G"] - 2.0) < 1e-12, "a missing value in one member does not erase the gene (mean of the others)")
empty = mc.concordance(coll[["GENEX"]], coll[["GENEY"]], {}, {}, mc.CONFIG_DEFAULTS)
check(empty["n_common"] == 0 and empty["note"] and empty["rows"] == [], "no gene in common says so, with no rows")
check(res["concordance"]["_matrices"][0].shape[1] == 300 and res["concordance"]["transcript_info"]["n_genes"] == 300, "the same-gene test uses ALL genes, not the top-variance subset")
tv2 = mc.run_correlation({"transcriptomics": T, "proteomics": Pm}, gene_of, {"top_variance": 20})
check(tv2["concordance"]["n_common"] == 250, "SENSITIVITY: lowering the top-variance cap does not shrink the same-gene analysis")

print("=== clustering and layout ===")
def naive_upgma(X):
    X = np.nan_to_num(np.asarray(X, float))
    n = len(X)
    D = [[math.sqrt(sum((X[i][k] - X[j][k]) ** 2 for k in range(X.shape[1]))) for j in range(n)] for i in range(n)]
    cl = {i: [i] for i in range(n)}
    nxt = n
    while len(cl) > 1:
        best = None
        keys = sorted(cl)
        for x in range(len(keys)):
            for y in range(x + 1, len(keys)):
                a_, b_ = keys[x], keys[y]
                d = sum(D[i][j] for i in cl[a_] for j in cl[b_]) / (len(cl[a_]) * len(cl[b_]))
                if best is None or d < best[0] - 1e-12:
                    best = (d, a_, b_)
        _, a_, b_ = best
        cl[nxt] = cl.pop(a_) + cl.pop(b_)
        nxt += 1
    return list(cl.values())[0]
rcl = np.random.RandomState(8)
blob1 = rcl.randn(5, 6) * 0.2 + 3
blob2 = rcl.randn(5, 6) * 0.2 - 3
shuf = np.r_[blob1, blob2][rcl.permutation(10)]
check(mc.cluster_order(shuf) == naive_upgma(shuf), "the clustering order equals a naive average-linkage written independently")
order = mc.cluster_order(shuf)
labels_ = [int(v.mean() > 0) for v in shuf]
seq = [labels_[i] for i in order]
check(sum(1 for i in range(9) if seq[i] != seq[i + 1]) == 1, f"SENSITIVITY: two planted groups end up contiguous ({seq})")
check(sorted(mc.cluster_order(shuf)) == list(range(10)) and mc.cluster_order(np.zeros((2, 3))) == [0, 1] and mc.cluster_order(np.zeros((1, 3))) == [0], "the order is a permutation; 1 and 2 rows are trivial")
cols = [[f"a{i}" for i in range(6)], [f"b{i}" for i in range(6)]]
edges = [(f"a{i}", f"b{5 - i}") for i in range(6)]
naive = {n: (c, r) for c, nodes in enumerate(cols) for r, n in enumerate(nodes)}
lay = mc.layered_layout(cols, edges)
check(mc.count_crossings(naive, edges) == 15 and mc.count_crossings(lay, edges) == 0,
      f"SENSITIVITY: the layered layout removes the crossings of a reversed matching (15 -> {mc.count_crossings(lay, edges)})")
check(mc.layered_layout(cols, edges) == lay and {c for c, _ in lay.values()} == {0, 1} and sorted(r for c, r in lay.values() if c == 0) == list(range(6)),
      "the layout is deterministic, keeps each node in its omics' column and uses distinct rows")
check(mc.count_crossings({"u": (0, 0), "v": (1, 1), "w": (0, 1), "x": (1, 0)}, [("u", "v"), ("w", "x")]) == 1, "crossings counted on a minimal example")
check(mc.layered_layout([["a"], []], []) == {"a": (0, 0)}, "a graph with no edges still places its nodes")

print("=== scatter vectors and labels ===")
pv = mc.pair_vectors(T, "t0", Mm, "m0")
check(pv["samples"] == idx and pv["x"] == [float(v) for v in T["t0"]] and len(pv["y"]) == n, "the scatter vectors are the two features' values per sample")
Tn = T.copy(); Tn.iloc[2, 0] = np.nan
check(mc.pair_vectors(Tn, "t0", Mm, "m0")["x"][2] is None and raises(mc.pair_vectors, T, "nope", Mm, "m0")[0] == "KeyError", "a missing value is None; an unknown feature is an error")
lab = mc.feature_labels(["a", "b", "c"], {"a": "TP53", "b": "TP53", "c": ""})
check(lab == {"a": "TP53", "b": "TP53 (2)", "c": "c"}, "labels use the gene symbol where there is one, are unique, and fall back to the id")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
