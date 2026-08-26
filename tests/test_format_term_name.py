# -*- coding: utf-8 -*-
"""Unit + integration tests for Genomics_module.format_term_name.

User-reported 2026-08-26: pathway/term names render lowercase-first
("cytoplasmic translation") across every module's tables and figures. Root
cause: GO's OWN naming convention is lowercase-first (measured live: 88.4%
of 21,028 real cached GO terms), and every enrichment engine passed each
database's published name straight through with no formatting at all.

A naive capitalize() is NOT safe here - 787 real names across GO, Reactome,
InterPro, UniProt keywords, Complex Portal and DISEASES start with a
scientifically load-bearing lowercase token (mRNA, cAMP, pH, eIF3, hnRNP,
bZIP, ...) that capitalizing would corrupt into a DIFFERENT, wrong token
(cAMP -> CAMP is a real, different human gene; pH -> PH is a real, different
protein domain). This file proves the fix both in isolation (the pure
function, against the REAL cached corpus, not just a fixture) and wired in
(the two shared per-namespace engines every Genomics/Transcriptomics/
Proteomics enrichment analysis routes through).
"""
import io
import json
import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as genomics_ora  # noqa: E402
import Transcriptomics_module as transcriptomics_engine  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


fmt = genomics_ora.format_term_name

# =============================================================================
# Basic behaviour
# =============================================================================
print("=== basic capitalization ===")
check(fmt("cytoplasmic translation") == "Cytoplasmic translation",
      f"the exact user-reported example is fixed (got {fmt('cytoplasmic translation')!r})")
check(fmt("Already Capitalized") == "Already Capitalized", "already-capitalized names are unchanged")
check(fmt("cell adhesion") == "Cell adhesion", "plain lowercase name capitalizes")
check(fmt("") == "", "empty string returns empty string, not an error")
check(fmt(None) is None, "None returns None, not an error")
check(fmt("a") == "A", "single-character name capitalizes")

# =============================================================================
# Hazard tokens: names capitalizing would corrupt into a DIFFERENT, real
# scientific token. Each of these is genuinely load-bearing, not synthetic -
# see the docstring above for cAMP/CAMP and pH/PH specifically.
# =============================================================================
print("=== lowercase-significant tokens are preserved verbatim ===")
hazard_examples = [
    "mRNA catabolic process", "tRNA processing", "rRNA processing",
    "snRNA processing", "lncRNA localization", "miRNA binding",
    "siRNA binding", "cAMP biosynthetic process", "cGMP binding",
    "dATP metabolic process", "dTDP-glucose 4,6-dehydratase",
    "pH reduction", "p53 binding", "p38MAPK cascade",
    "eIF3 complex", "eIF2B Complex Subunit Epsilon",
    "hnRNP A1, RNA recognition motif 1", "bZIP transcription factor complex",
    "bHLH transcription factor binding", "eRF1 domain 2",
]
for name in hazard_examples:
    out = fmt(name)
    check(out == name, f"{name!r} preserved verbatim (got {out!r})")

# =============================================================================
# MSigDB SCREAMING_SNAKE_CASE: underscores -> spaces, case fully preserved
# (TNFA/NFKB are real gene symbols a case change would corrupt).
# =============================================================================
print("=== MSigDB underscore handling ===")
check(fmt("HALLMARK_TNFA_SIGNALING_VIA_NFKB") == "HALLMARK TNFA SIGNALING VIA NFKB",
      f"underscores become spaces, caps untouched (got {fmt('HALLMARK_TNFA_SIGNALING_VIA_NFKB')!r})")
check(fmt("GOBP_MRNA_SPLICING") == "GOBP MRNA SPLICING", "another MSigDB collection prefix behaves the same")
check(fmt("KEGG") == "KEGG", "an all-caps name with no underscore is left alone (no-op replace)")
# A mixed-case name that happens to contain an underscore must NOT hit the
# MSigDB branch - it should fall through to the normal lowercase-guard path.
check(fmt("biological_process") == "Biological_process",
      f"a lowercase name with an underscore (GO's own root term) is NOT treated as "
      f"MSigDB-style - only capitalized, underscore left alone (got {fmt('biological_process')!r})")

# =============================================================================
# Idempotence: running the formatter twice must be a no-op the second time -
# every call site applies it exactly once, but this guards against a future
# double-application anywhere not silently mangling a name further.
# =============================================================================
print("=== idempotence ===")
for name in ["cell adhesion", "mRNA catabolic process", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Already Capitalized", ""]:
    once, twice = fmt(name), fmt(fmt(name))
    check(once == twice, f"fmt is idempotent for {name!r} (got {once!r} then {twice!r})")

# =============================================================================
# SENSITIVITY: recompute over the REAL cached universes this app actually
# uses, not a hand-picked fixture - proves the fix's real-world impact and
# pins it as a regression guard. If a future edit to the regex silently
# changes behaviour on the real corpus, this catches it.
# =============================================================================
print("=== sensitivity: full real-corpus verification ===")


def _names_from(path, key=None):
    d = json.load(io.open(path, encoding="utf-8"))
    out = []
    if key:
        for info in d.get(key, {}).values():
            out.append(info.get("name", ""))
    else:
        for ns in ("biological_process", "molecular_function", "cellular_component"):
            for info in d.get(ns, {}).values():
                out.append(info.get("name", ""))
    return [n for n in out if n]


go_path = os.path.join(_ROOT, "cache", "ora_universe", "GO_v2_9606.json")
if os.path.isfile(go_path):
    go_names = _names_from(go_path)
    lower_start = [n for n in go_names if n[:1].islower()]
    transformed = [fmt(n) for n in lower_start]
    n_capitalized = sum(1 for old, new in zip(lower_start, transformed) if old != new)
    n_preserved = len(lower_start) - n_capitalized
    check(len(lower_start) > 15000,
          f"sanity: the real cached GO universe still has a large lowercase-start majority "
          f"(got {len(lower_start)} of {len(go_names)} total - if this ever drops near 0 the "
          f"cache was probably rebuilt with different content, re-verify the numbers below)")
    # Measured live 2026-08-26: 18,580 lowercase-start GO terms, 352 hazard-
    # preserved (18,228 capitalized). Pin these exact counts so any change to
    # the regex is caught, not just "still positive".
    check(n_capitalized == 18228,
          f"exactly 18228 real GO terms get capitalized (got {n_capitalized}) - "
          f"if this moved, the hazard regex's coverage changed")
    check(n_preserved == 352,
          f"exactly 352 real GO terms are hazard-preserved (got {n_preserved})")
    # Zero false positives: every preserved name must genuinely start with a
    # lowercase-then-uppercase/digit run, not something the regex over-matched.
    preserved = [n for n, t in zip(lower_start, transformed) if n == t]
    false_positive_re = re.compile(r"^[a-z]{1,3}[A-Z0-9]")
    non_matching = [n for n in preserved if not false_positive_re.match(n)]
    check(len(non_matching) == 0,
          f"every preserved name genuinely matches the hazard shape - zero false positives "
          f"(got {len(non_matching)} unexplained: {non_matching[:5]})")
else:
    print("  (GO cache not present - skipping full-corpus sensitivity check, ran the fixture-based checks above instead)")

# =============================================================================
# SENSITIVITY: prove a naive single-character capitalize() - the realistic
# accidental implementation this fix replaced - would genuinely corrupt real
# names, so the guard logic above is not solving an imaginary problem.
# =============================================================================
print("=== SENSITIVITY: the naive approach this replaces really is wrong ===")


def _naive_capitalize(s):
    return s[0].upper() + s[1:] if s else s


naive_mrna = _naive_capitalize("mRNA catabolic process")
check(naive_mrna == "MRNA catabolic process",
      f"confirmed: naive capitalize turns 'mRNA' into the wrong token 'MRNA' (got {naive_mrna!r})")
check(fmt("mRNA catabolic process") != naive_mrna,
      "the real formatter's output differs from the naive approach on this exact input - "
      "proving the guard is load-bearing, not decorative")
naive_ph = _naive_capitalize("pH reduction")
check(naive_ph == "PH reduction",
      f"confirmed: naive capitalize turns chemistry 'pH' into 'PH', the pleckstrin-homology "
      f"domain - a different, real, wrong entity (got {naive_ph!r})")

# =============================================================================
# INTEGRATION: the two shared per-namespace engines every Genomics ORA,
# Genomics Burden Test, Transcriptomics ORA, and (by reuse) Proteomics
# ORA/Quantitative Enrichment analysis routes through. Prove format_term_name
# is actually WIRED IN, not just correct in isolation.
# =============================================================================
print("=== INTEGRATION: Genomics_module._run_fisher_for_namespace (ORA) ===")
term_map = {
    "T1": {"name": "cell adhesion", "genes": {"A", "B", "C"}},
    "T2": {"name": "mRNA catabolic process", "genes": {"D", "E", "F"}},
}
universe = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"}
sig_set = {"A", "B", "D", "E"}
block = genomics_ora._run_fisher_for_namespace(term_map, sig_set, universe, 0.05, True, 1, None)
names_out = dict(zip(block["table"]["Term ID"], block["table"]["Term Name"]))
check(names_out.get("T1") == "Cell adhesion", f"ORA: normal name capitalized in real output (got {names_out.get('T1')!r})")
check(names_out.get("T2") == "mRNA catabolic process",
      f"ORA: hazard name preserved in real output, not corrupted to 'MRNA...' (got {names_out.get('T2')!r})")

print("=== INTEGRATION: Genomics_module._run_wilcoxon_for_namespace (Burden Test / Proteomics Quant) ===")
score_by_gene = {g: (1.0 if g in {"A", "B", "C"} else -1.0) for g in universe}
wblock = genomics_ora._run_wilcoxon_for_namespace(term_map, score_by_gene, universe, True, 1, None, alternative="two-sided")
wnames_out = dict(zip(wblock["table"]["Term ID"], wblock["table"]["Term Name"]))
check(wnames_out.get("T1") == "Cell adhesion", f"Burden/Quant: normal name capitalized (got {wnames_out.get('T1')!r})")
check(wnames_out.get("T2") == "mRNA catabolic process",
      f"Burden/Quant: hazard name preserved (got {wnames_out.get('T2')!r})")

print("=== INTEGRATION: Transcriptomics_module._run_ora_for_sets (Transcriptomics ORA) ===")
t_term_map = {
    "TT1": {"name": "cell adhesion", "genes": {"GENE_A", "GENE_B", "GENE_C"}},
    "TT2": {"name": "tRNA processing", "genes": {"GENE_D", "GENE_E", "GENE_F"}},
}
t_universe = {"GENE_" + c for c in "ABCDEFGHIJ"}
t_degs = {"GENE_A", "GENE_B", "GENE_D", "GENE_E"}
t_block = transcriptomics_engine._run_ora_for_sets(t_term_map, t_universe, t_degs, 1, None, None)
t_names_out = dict(zip(t_block["table"]["Term ID"], t_block["table"]["Term Name"]))
check(t_names_out.get("TT1") == "Cell adhesion", f"Transcriptomics ORA: normal name capitalized (got {t_names_out.get('TT1')!r})")
check(t_names_out.get("TT2") == "tRNA processing",
      f"Transcriptomics ORA: hazard name preserved (got {t_names_out.get('TT2')!r})")

print("=== INTEGRATION (source-level): Transcriptomics_module.run_gsea wires format_term_name too ===")
# GSEA's row-builder needs a full ranked-list + permutation run to reach
# offline (no minimal synthetic entry point exists, unlike the ORA/Burden
# engines above) - checked at the source level instead: confirm the exact
# line format_term_name wraps o["name"], so a future edit cannot silently
# drop the wrap without this failing.
gsea_src = io.open(os.path.join(_ROOT, "Transcriptomics_module.py"), encoding="utf-8").read()
check('"Term Name": format_term_name(o["name"])' in gsea_src,
      "run_gsea's row-builder wraps the term name in format_term_name (source-level check)")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    sys.exit(1)
