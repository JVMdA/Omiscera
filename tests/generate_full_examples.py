"""One-off generator: build 3000-gene Ensembl-ID and Entrez-ID versions of
Data3/example_transcriptomics_de_full.csv, so every identifier type's
"Example file" button produces a dataset large enough to actually run TF/
Regulon activity inference (and GSEA) against - not just gene_symbol.

Measured problem this fixes: the ensembl/entrez examples are 18 genes each,
against a default minimum regulon size of 5, so essentially no TF clears the
bar and the activity table comes back EMPTY.

Reuses the app's own retriever (resolve_de_table) rather than an ad-hoc
mapping, so the IDs written here are exactly the ones the app would resolve
back. Reports survivor counts - a silently shrunken example would be its own
quiet defect.

Transcript is deliberately NOT regenerated: that example exists to
demonstrate transcript COLLAPSING (a retriever behaviour), and enlarging it
would break what it is for.
"""
import os
import sys

import pandas as pd

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_module_retriever as retr

DATA3 = os.path.join(_ROOT, "Data3")
SRC = DATA3 + r"\example_transcriptomics_de_full.csv"

print(f"Reading {SRC} ...")
with open(SRC, "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
n_input = len(de_df)
print(f"  input rows: {n_input}")

print("Resolving via the app's own retriever (mygene.info) ...")
resolved_df, meta = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)
print(f"  resolved rows: {len(resolved_df)}")
print(f"  mapped={meta.get('n_mapped')} unmapped={meta.get('n_unmapped')}")


def _write(target_col, out_name, header_label):
    df = resolved_df.copy()
    df[target_col] = df[target_col].astype(str).str.strip()
    # Keep only rows whose target ID space actually resolved, and drop
    # duplicates on that ID (a gene must appear once, exactly as the
    # retriever's own collapsing guarantees for the symbol example).
    df = df[(df[target_col] != "") & (df[target_col].str.lower() != "nan")]
    before_dedup = len(df)
    df = df.drop_duplicates(subset=[target_col], keep="first")
    out = pd.DataFrame({
        header_label: df[target_col],
        "log2FoldChange": df["log2FC"],
        "pvalue": df["pvalue"],
        "padj": df["padj"],
    })
    path = DATA3 + "\\" + out_name
    out.to_csv(path, index=False)
    pct = 100.0 * len(out) / n_input
    print(f"  {out_name}: {len(out)} rows written "
          f"({pct:.1f}% of the {n_input} input genes survived mapping; "
          f"{before_dedup - len(df)} duplicate IDs collapsed)")
    return len(out)


print("Writing new example files ...")
n_ens = _write("ensembl_id", "example_transcriptomics_ensembl_full.csv", "Ensembl Gene ID")
n_ent = _write("entrez_id", "example_transcriptomics_entrez_full.csv", "Entrez Gene ID")

print("\nDone.")
print(f"  Ensembl example: {n_ens} genes (was 18)")
print(f"  Entrez example:  {n_ent} genes (was 18)")
print("  Transcript example: intentionally left at 7 genes (demonstrates collapsing).")
