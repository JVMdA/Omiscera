"""Guards the Transcriptomics "Example file" downloads.

Two defects this exists to prevent recurring:
  1. API.py decides WHAT is served; static/app.js decides what the browser
     SAVES it as (downloadFile sets a.download). They had drifted:
     gene_symbol served the 3000-gene example_transcriptomics_de_full.csv but
     saved it as "example_transcriptomics_symbol.csv" - the name of a real,
     different 20-gene file. Nothing caught it because nothing compared the
     two maps.
  2. The ensembl/entrez examples were 18 genes each, far too small for
     TF/Regulon activity inference (min regulon size 5) or GSEA (min term
     size 15), so "Example file -> Resolve -> Run" returned an EMPTY table.

No server needed - parses both source files and inspects Data3 directly.
"""
import os
import re
import sys

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import API

BASE = _ROOT
DATA3 = os.path.join(BASE, "Data3")
APP_JS = os.path.join(BASE, "static", "app.js")

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def parse_js_example_map(path):
    """Extract TRANSCRIPTOMICS_EXAMPLE_FILENAMES from app.js."""
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    m = re.search(r"const TRANSCRIPTOMICS_EXAMPLE_FILENAMES\s*=\s*\{(.*?)\};", text, re.S)
    if not m:
        return None
    body = m.group(1)
    out = {}
    for key, val in re.findall(r"(\w+)\s*:\s*[\"']([^\"']+)[\"']", body):
        out[key] = val
    return out


print("=== the two example maps must agree ===")
js_map = parse_js_example_map(APP_JS)
py_map = API.TRANSCRIPTOMICS_EXAMPLE_BY_IDENTIFIER_TYPE
check(js_map is not None, "TRANSCRIPTOMICS_EXAMPLE_FILENAMES parsed out of app.js")
if js_map is not None:
    check(set(js_map.keys()) == set(py_map.keys()),
          f"both maps cover the same identifier types (js={sorted(js_map)}, py={sorted(py_map)})")
    for key in sorted(set(js_map) & set(py_map)):
        check(js_map[key] == py_map[key],
              f"{key}: served '{py_map[key]}' == saved-as '{js_map[key]}'")

print("=== every served example actually exists ===")
for key, fname in sorted(py_map.items()):
    path = os.path.join(DATA3, fname)
    check(os.path.isfile(path), f"{key}: {fname} exists in Data3")

print("=== examples are large enough for the downstream analyses ===")
# TF/Regulon needs >= REGULON min regulon size targets per TF present, and
# GSEA needs >= 15 genes per term - both are hopeless on ~18 genes. The three
# non-transcript examples must therefore be substantial.
MIN_USEFUL_ROWS = 1000
for key in ("gene_symbol", "ensembl_gene", "entrez"):
    fname = py_map[key]
    path = os.path.join(DATA3, fname)
    if not os.path.isfile(path):
        check(False, f"{key}: file missing, cannot count rows")
        continue
    with open(path, "r", encoding="utf-8") as fh:
        n_rows = sum(1 for _ in fh) - 1   # minus header
    check(n_rows >= MIN_USEFUL_ROWS,
          f"{key}: {fname} has {n_rows} genes (>= {MIN_USEFUL_ROWS}, enough for TF/Regulon + GSEA)")

# The transcript example is DELIBERATELY small - it demonstrates transcript
# collapsing, a retriever behaviour. Assert it stays small ON PURPOSE, so
# that intent is recorded rather than looking like the same oversight.
transcript_path = os.path.join(DATA3, py_map["ensembl_transcript"])
with open(transcript_path, "r", encoding="utf-8") as fh:
    n_transcript = sum(1 for _ in fh) - 1
check(n_transcript < 50,
      f"ensembl_transcript is intentionally small ({n_transcript} rows) - it exists to demonstrate "
      f"transcript collapsing, not to feed downstream analyses")

print("=== example headers match the identifier type they claim ===")
EXPECTED_HEADER = {
    "gene_symbol": "Gene Symbol",
    "ensembl_gene": "Ensembl Gene ID",
    "ensembl_transcript": "Ensembl Transcript ID",
    "entrez": "Entrez Gene ID",
}
for key, fname in sorted(py_map.items()):
    path = os.path.join(DATA3, fname)
    if not os.path.isfile(path):
        continue
    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline().strip()
    first_col = header.split(",")[0].strip().lstrip("\ufeff")
    check(first_col == EXPECTED_HEADER[key],
          f"{key}: first column is '{first_col}' (expected '{EXPECTED_HEADER[key]}')")
    check("log2FoldChange" in header and "pvalue" in header,
          f"{key}: header carries log2FoldChange + pvalue")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
