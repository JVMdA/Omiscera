# -*- coding: utf-8 -*-
"""Unit tests for Metabolomics_module's concentration-table pipeline
(parse_concentration_table, compute_compound_statistics,
build_scores_from_stats), offline.

THE CORRECTION THIS FILE VERIFIES: an earlier version of this module wrongly
asked users to upload a pre-computed log2FC/pvalue/padj per compound - the
RNA-seq/proteomics differential-expression convention, which is not what a
metabolomics researcher actually has. MetaboAnalyst's own Quantitative
Enrichment Analysis (QEA) takes a raw concentration table and computes
statistics itself; this file is the same design, and these tests are what
would have caught the original mistake had they existed first.
"""
import io
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Metabolomics_module as mm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def _csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


# =============================================================================
print("=== parse_concentration_table: standard orientation ===")
# =============================================================================
rng = np.random.default_rng(0)
rows = []
for i in range(5):
    rows.append({"Sample": f"C{i}", "Group": "Control",
                "Citrate": float(rng.normal(10, 1)), "Succinate": float(rng.normal(5, 0.5))})
for i in range(5):
    rows.append({"Sample": f"T{i}", "Group": "Treated",
                "Citrate": float(rng.normal(40, 4)), "Succinate": float(rng.normal(5, 0.5))})
standard_df = pd.DataFrame(rows)

values_df, groups, meta = mm.parse_concentration_table(_csv_bytes(standard_df))
check(meta["orientation"] == "standard", f"correctly detected as standard orientation (got {meta['orientation']!r})")
check(meta["group_column"] == "Group", f"group column found by name (got {meta['group_column']!r})")
check(sorted(meta["group_levels"]) == ["Control", "Treated"],
      f"both group levels detected (got {meta['group_levels']})")
check(values_df.shape == (10, 2), f"values matrix is samples x compounds (got {values_df.shape})")
check(list(values_df.columns) == ["Citrate", "Succinate"],
      f"compound columns preserved in order (got {list(values_df.columns)})")
check(groups.tolist().count("Control") == 5 and groups.tolist().count("Treated") == 5,
      "group assignment matches the uploaded labels exactly")

# =============================================================================
print("=== parse_concentration_table: transposed orientation (auto-detected) ===")
# =============================================================================
# Fixture: compounds in rows, samples in columns, with the group labels carried as
# a pseudo-compound row named "Group" (which becomes the group column once transposed).
wide_rows = [
    ["Group"] + ["Control"] * 5 + ["Treated"] * 5,
    ["Citrate"] + [10.1, 9.8, 10.5, 9.6, 10.2, 40.1, 38.5, 41.2, 39.8, 40.5],
    ["Succinate"] + [5.0, 5.1, 4.9, 5.2, 5.0, 5.1, 4.8, 5.3, 5.0, 4.9],
]
wide_df = pd.DataFrame(wide_rows, columns=["Compound"] + [f"S{i}" for i in range(10)])
values_wide, groups_wide, meta_wide = mm.parse_concentration_table(_csv_bytes(wide_df))
check(meta_wide["orientation"] == "transposed",
      f"a compounds-in-rows file is auto-detected as transposed (got {meta_wide['orientation']!r})")
check(meta_wide["group_column"] == "Group", "the 'Group' pseudo-compound row becomes the group column after transposing")
check(values_wide.shape == (10, 2), f"values matrix ends up samples x compounds either way (got {values_wide.shape})")
check(groups_wide.tolist().count("Control") == 5 and groups_wide.tolist().count("Treated") == 5,
      "group labels survive the transpose correctly")

# =============================================================================
print("=== parse_concentration_table: group column detection fallback ===")
# =============================================================================
unnamed = standard_df.rename(columns={"Group": "Cohort"})
values_u, groups_u, meta_u = mm.parse_concentration_table(_csv_bytes(unnamed))
check(meta_u["group_column"] == "Cohort",
      f"an unrecognised-but-plausible label column (2-10 levels, >=2 per level) is still found "
      f"(got {meta_u['group_column']!r})")

# =============================================================================
print("=== parse_concentration_table: honest failures ===")
# =============================================================================
too_few_cols = pd.DataFrame({"Sample": ["A", "B"], "Group": ["X", "Y"]})
try:
    mm.parse_concentration_table(_csv_bytes(too_few_cols))
    check(False, "a table with no numeric compound column should raise")
except ValueError as exc:
    check("at least" in str(exc).lower(), f"clear error naming the missing column count (got: {exc})")

one_group = standard_df.copy()
one_group["Group"] = "OnlyOne"
try:
    mm.parse_concentration_table(_csv_bytes(one_group))
    check(False, "a table with only ONE group should raise - nothing to compare")
except ValueError as exc:
    check("two" in str(exc).lower() or "group" in str(exc).lower(),
          f"error explains at least two groups are needed (got: {exc})")

no_group_col = pd.DataFrame({
    "Sample": [f"S{i}" for i in range(6)],
    "Compound_A": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    "Compound_B": [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
})
try:
    mm.parse_concentration_table(_csv_bytes(no_group_col))
    check(False, "a table with no plausible group column at all should raise")
except ValueError as exc:
    check("group" in str(exc).lower() or "class" in str(exc).lower(),
          f"error explains a group/class column is needed (got: {exc})")

# =============================================================================
print("=== compute_compound_statistics: recovers a known planted difference ===")
# =============================================================================
values_df, groups, _ = mm.parse_concentration_table(_csv_bytes(standard_df))
stats_df, stats_meta = mm.compute_compound_statistics(values_df, groups, "Treated", "Control")
check(len(stats_df) == 2, f"both compounds tested (got {len(stats_df)})")
citrate_row = stats_df[stats_df["Name"] == "Citrate"].iloc[0]
succinate_row = stats_df[stats_df["Name"] == "Succinate"].iloc[0]
check(citrate_row["log2fc"] > 1.5,
      f"Citrate's planted ~4x shift (Treated >> Control) recovers a strongly positive log2fc "
      f"(got {citrate_row['log2fc']:.3f})")
check(citrate_row["pvalue"] < 0.001,
      f"Citrate's clean, large, well-separated shift is genuinely significant "
      f"(got p={citrate_row['pvalue']:.3g})")
check(succinate_row["pvalue"] > 0.05,
      f"Succinate has NO planted difference and correctly does not reach significance "
      f"(got p={succinate_row['pvalue']:.3g})")
check("padj" in stats_df.columns and stats_df["padj"].notna().all(),
      "BH-adjusted p-values are computed for every tested compound")
check(stats_meta["group_a"] == "Treated" and stats_meta["group_b"] == "Control",
      "meta records which group is which side of the contrast")
check(stats_meta["log_transform"] is True and stats_meta["log_offset"] > 0,
      f"log-transform is on by default with a positive, reported offset "
      f"(got offset={stats_meta['log_offset']})")

# =============================================================================
print("=== compute_compound_statistics: min_per_group drops honestly, no imputation ===")
# =============================================================================
sparse = standard_df.copy()
sparse.loc[0:2, "Citrate"] = np.nan  # only 2 of 5 Control values remain for Citrate
sparse_values, sparse_groups, _ = mm.parse_concentration_table(_csv_bytes(sparse))
sparse_stats, sparse_meta = mm.compute_compound_statistics(
    sparse_values, sparse_groups, "Treated", "Control", min_per_group=3)
check("Citrate" not in set(sparse_stats["Name"]),
      "a compound with fewer than min_per_group finite values in a group is DROPPED, "
      "never imputed or silently kept with a degenerate test")
check(sparse_meta["n_dropped_low_n"] == 1,
      f"the drop is counted, not silent (got n_dropped_low_n={sparse_meta['n_dropped_low_n']})")
check("Succinate" in set(sparse_stats["Name"]),
      "an unaffected compound in the same upload is unaffected by another compound's missingness")

# =============================================================================
print("=== compute_compound_statistics: too few samples in a group raises clearly ===")
# =============================================================================
# One Control and one Treated sample: build the matrix directly (the parser itself
# already rejects a group with fewer than two samples, so go around it).
one_each = pd.DataFrame({"Citrate": [1.0, 2.0], "Succinate": [3.0, 4.0]}, index=["C0", "T0"])
one_each_groups = pd.Series(["Control", "Treated"], index=one_each.index)
try:
    mm.compute_compound_statistics(one_each, one_each_groups, "Treated", "Control")
    check(False, "fewer than 2 samples in a group should raise")
except ValueError as exc:
    check("2 samples" in str(exc) or "at least 2" in str(exc).lower(),
          f"error names the minimum-samples requirement clearly (got: {exc})")

# =============================================================================
print("=== build_scores_from_stats: joins by Name, mirrors the standing Wilcoxon contract ===")
# =============================================================================
matched = pd.DataFrame([
    {"Name": "Citrate", "compound_key": "KEY_CITRATE", "match_tier": "inchikey", "class_level_only": False},
    {"Name": "Succinate", "compound_key": "KEY_SUCCINATE", "match_tier": "inchikey", "class_level_only": False},
    {"Name": "Unmeasured Compound", "compound_key": "", "match_tier": "", "class_level_only": False},
], columns=mm.MATCHED_COLUMNS)

scores, n_dropped = mm.build_scores_from_stats(matched, stats_df, metric="log2fc")
check(scores.get("KEY_CITRATE", 0) > 1.5,
      f"Citrate's score under the default log2fc metric matches its computed statistic "
      f"(got {scores.get('KEY_CITRATE')})")
check("KEY_SUCCINATE" in scores, "every tested compound with a matched key gets a score")
check(n_dropped == 0, f"nothing unexpectedly dropped in this fixture (got n_dropped={n_dropped})")

scores_t, _ = mm.build_scores_from_stats(matched, stats_df, metric="t_stat")
check(abs(scores_t["KEY_CITRATE"] - citrate_row["t_stat"]) < 1e-9,
      "metric='t_stat' selects Welch's t rather than log2fc")

bad_metric_scores, _ = mm.build_scores_from_stats(matched, stats_df, metric="not_a_real_metric")
check(bad_metric_scores == scores,
      "an unrecognised metric falls back to QUANT_METRIC_DEFAULT rather than raising or "
      "silently returning nothing")

# =============================================================================
print("=== explicit overrides are honoured or rejected, never silently replaced ===")
# =============================================================================
_, _, meta_cohort = mm.parse_concentration_table(_csv_bytes(unnamed))
check("Cohort" in meta_cohort["candidate_group_columns"],
      f"a group column that is not NAMED like one but looks like one is offered as a choice "
      f"(got {meta_cohort['candidate_group_columns']})")
try:
    mm.parse_concentration_table(_csv_bytes(standard_df), group_column="Nope")
    check(False, "an unknown group_column override should raise, not fall back to auto-detection")
except ValueError as exc:
    check("Nope" in str(exc), f"the error names the column that was not found (got: {exc})")
try:
    mm.parse_concentration_table(_csv_bytes(standard_df), group_column="Sample")
    check(False, "the sample-ID column cannot be the group column")
except ValueError:
    check(True, "the sample-ID column is rejected as a group column")
_, groups_forced, meta_forced = mm.parse_concentration_table(_csv_bytes(unnamed), group_column="Cohort")
check(meta_forced["group_column"] == "Cohort" and sorted(meta_forced["group_levels"]) == ["Control", "Treated"],
      "a valid explicit group_column is used as given")

print("=== a compound listed TWICE (two features of one metabolite) is kept, not a crash ===")
# Found with a real table: one metabolite often appears as two features (two
# peaks / adducts). Selecting by label returned every match per occurrence, and
# a repeated label made the group-column check receive a DataFrame - a crash.
_S = [f"E{e}S{i}" for e in (1, 2) for i in range(1, 7)]
_rng = np.random.default_rng(5)


def _row(name, shift):
    v = _rng.normal(0, 1, 12)
    v[6:] += shift
    return f'"{name}",' + ",".join(f"{x:.4f}" for x in v)


_dup = "\n".join(['"",' + ",".join(_S), '"Label",' + ",".join(x[:2] for x in _S),
                   _row("Malic acid", 2), _row("Malic acid", 3), _row("Citric acid", 0)]).encode()
_dv, _dg, _dm = mm.parse_concentration_table(_dup)
check(list(_dv.columns) == ["Malic acid", "Malic acid", "Citric acid"] and _dv.shape == (12, 3),
      f"both Malic acid features are kept as separate columns, none multiplied (got {list(_dv.columns)})")
check(_dm["n_repeated_names"] == 1 and _dm["n_compounds"] == 3, "meta reports the repeated name")
_dst, _ = mm.compute_compound_statistics(_dv, _dg, "E1", "E2", log_transform=False)
_mal = _dst[_dst["Name"] == "Malic acid"]["log2fc"].tolist()
check(len(_dst) == 3 and len(_mal) == 2 and abs(_mal[0] - _mal[1]) > 0.5,
      f"each feature gets its OWN statistics, not the first one twice (got {_mal})")
_dmatched = pd.DataFrame({"Name": ["Malic acid", "Citric acid"], "compound_key": ["K1", "K2"]})
_dscores, _ = mm.build_scores_from_stats(_dmatched, _dst, "log2fc")
check(abs(_dscores["K1"] - sum(_mal) / 2) < 1e-9, "the two features are averaged into ONE score for the compound")

print("=== blank rows/columns and unnamed or all-empty compounds are ignored, and counted ===")
_g = [f"A{i}" for i in range(1, 7)] + [f"B{i}" for i in range(1, 7)]
_r2 = np.random.default_rng(9)
_f = lambda: ",".join(f"{x:.3f}" for x in _r2.normal(size=12))
# STANDARD layout: a blank column, an all-missing compound, trailing blank columns/rows
_std = "Sample,Group,Citrate,,Glucose,AllMissing,,\n" + "\n".join(
    f"{s},{s[0]},{_r2.normal():.3f},,{_r2.normal():.3f},,," for s in _g) + "\n,,,,,,,\n,,,,,,,\n"
_sv, _sg, _sm = mm.parse_concentration_table(_std.encode())
check(list(_sv.columns) == ["Citrate", "Glucose"] and _sv.shape == (12, 2),
      f"standard layout: only the real compounds remain (got {list(_sv.columns)}, shape {_sv.shape})")
check(_sm["n_compounds"] == 2 and sorted(_sm["group_levels"]) == ["A", "B"] and _sm["n_samples"] == 12,
      f"...blank trailing rows are not samples or a bogus group (levels {_sm['group_levels']}, n={_sm['n_samples']})")
check(_sm["n_empty_dropped"] == 1, f"the all-missing compound is counted as ignored (got {_sm['n_empty_dropped']})")
# TRANSPOSED layout: a named row of data with a BLANK name, an all-empty named row, blank rows
_tr = "\n".join(["Sample," + ",".join(_g), "Group," + ",".join(x[0] for x in _g),
                  "Citrate," + _f(), "," + _f(), "Glucose," + _f(),
                  "EmptyRow," + "," * 11, "," + "," * 11, "," + "," * 11])
_tv, _tg, _tm = mm.parse_concentration_table(_tr.encode())
check(list(_tv.columns) == ["Citrate", "Glucose"],
      f"transposed layout: the unnamed and the empty rows are ignored (got {list(_tv.columns)})")
check(_tm["n_unnamed_dropped"] == 1 and _tm["n_empty_dropped"] == 1,
      f"...and each is counted (unnamed {_tm['n_unnamed_dropped']}, empty {_tm['n_empty_dropped']})")
_te, _ = mm.compute_compound_statistics(_tv, _tg, "A", "B", log_transform=False)
check("nan" not in _te["Name"].astype(str).str.lower().tolist(), "no compound called 'nan' reaches the statistics")

print("=== negative values with log2 ticked are refused, never silently clipped ===")
_ni = [f"A{i}" for i in range(6)] + [f"B{i}" for i in range(6)]
_nv = pd.DataFrame({"X": np.r_[np.full(6, -1.0) + np.arange(6) * 0.05, np.full(6, 1.0) + np.arange(6) * 0.05],
                    "Y": np.linspace(-1, 1, 12)}, index=_ni)
_ng = pd.Series([i[0] for i in _ni], index=_ni)
try:
    mm.compute_compound_statistics(_nv, _ng, "A", "B", log_transform=True)
    check(False, "log2 on negative values must raise")
except ValueError as exc:
    check("negative" in str(exc) and "Untick" in str(exc),
          f"it raises a ValueError that names the problem AND the fix (got {str(exc)[:80]!r}...)")
_nst, _ = mm.compute_compound_statistics(_nv, _ng, "A", "B", log_transform=False)
check(abs(float(_nst.loc[_nst.Name == "X", "log2fc"].iloc[0]) - (-2.0)) < 1e-9,
      "with log2 off the mean difference is exact (-2.0), not the clipped value the old code returned (-2.8)")
_pos = _nv.abs() + 1
check(len(mm.compute_compound_statistics(_pos, _ng, "A", "B", log_transform=True)[0]) == 2,
      "ordinary positive data with log2 ticked is unaffected")

print("=== layout is detected by STRUCTURE, so the group label can be called anything ===")
_o = np.random.default_rng(11)
_sm = [f"S{i}" for i in range(1, 13)]
_gl = ["Ctrl"] * 6 + ["Trt"] * 6
_cmp = {"Citrate": _o.normal(10, 1, 12), "Malate": _o.normal(8, 1, 12), "Glucose": _o.normal(5, 1, 12)}
def _std_table(group_header):
    lines = ["Sample," + group_header + "," + ",".join(_cmp)]
    for i, s in enumerate(_sm):
        lines.append(f"{s},{_gl[i]}," + ",".join(f"{_cmp[k][i]:.3f}" for k in _cmp))
    return "\n".join(lines).encode()
def _tr_table(group_header):
    lines = ["Sample," + ",".join(_sm), group_header + "," + ",".join(_gl)]
    for k in _cmp:
        lines.append(k + "," + ",".join(f"{x:.3f}" for x in _cmp[k]))
    return "\n".join(lines).encode()
for _hdr in ["Group", "Sample type", "Cohort", "Arm", "Study_arm"]:
    _v1, _g1, _m1 = mm.parse_concentration_table(_std_table(_hdr))
    check(_m1["orientation"] == "standard" and list(_v1.columns) == list(_cmp) and sorted(_m1["group_levels"]) == ["Ctrl", "Trt"],
          f"samples in rows, group column called {_hdr!r}: read as standard with 3 compounds")
    _v2, _g2, _m2 = mm.parse_concentration_table(_tr_table(_hdr))
    check(_m2["orientation"] == "transposed" and list(_v2.columns) == list(_cmp) and sorted(_m2["group_levels"]) == ["Ctrl", "Trt"],
          f"compounds in rows, group row called {_hdr!r}: read as transposed with 3 compounds (got {_m2['orientation']})")
_a, _, _ = mm.parse_concentration_table(_std_table("Cohort"))
_b, _, _ = mm.parse_concentration_table(_tr_table("Cohort"))
check(np.allclose(_a.values, _b.values), "the two layouts of the same data give identical numbers")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
