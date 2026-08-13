TRANSCRIPTOMICS — RETRIEVER
=======================================

WHAT THIS STEP IS FOR
--------------------------
Before any Transcriptomics analysis can run, your differential-expression
(DE) results have to be turned into one clean, consistent table: every
gene identified the same way, carrying its expression signal (log2
fold-change and p-value), with nothing thresholded away and nothing
silently dropped. This is the shared first step for every Transcriptomics
analysis — all of them read the same resolved table this step produces.

The retriever performs no statistics itself. It doesn't threshold your
data, rank it, or compare it against anything — it only checks your file,
matches gene identifiers, and hands back your full, unmodified expression
signal for the analysis stage to use.

USER INPUT
-------------
A differential-expression results table (CSV/TSV) — the typical output of
a DESeq2, edgeR, or limma run:

  - One identifier column (required) — pick which type you're uploading:
      - Gene Symbol
      - Ensembl Gene ID
      - Ensembl Transcript ID
      - Entrez Gene ID
  - log2FoldChange (required) — your expression signal.
  - pvalue (required).
  - An adjusted p-value / FDR column is OPTIONAL — if missing, it's
    computed automatically (Benjamini-Hochberg) so every downstream
    analysis always has one to work with.

Column headers don't need to match exactly — common naming variants
(log2FC, logFC, LFC; padj, FDR, qvalue; and similar) are recognized
automatically, and the screen shows you which uploaded column was matched
to which required field.

An organism (species) selection is always required, regardless of
identifier type.

IDENTIFIER MAPPING
----------------------
Whatever identifier type you upload, every gene is mapped to its official
gene symbol (plus Entrez and Ensembl gene IDs, when available), so the
table can be matched against any downstream gene-set database, whatever
identifier space that database uses.

This works for ANY organism your selected species corresponds to — not
just human or mouse.

Each mapped gene gets a status:
  - MAPPED — resolved directly.
  - ALIAS_MAPPED — resolved via an older/alternative gene name rather than
    its current official symbol (e.g. "HER2" resolves to "ERBB2"). Purely
    informational — treated identically to a direct match afterward.
  - TRANSCRIPT_COLLAPSED — see below.
  - UNMAPPED — no match found. Unmapped rows are always kept in the
    output, never silently dropped, so you can see exactly what wasn't
    resolved.

TRANSCRIPT-LEVEL INPUT
---------------------------
If you upload Ensembl Transcript IDs, each transcript is first mapped to
its parent gene. Since a gene can have several transcripts in your table,
all transcripts mapping to the same gene are collapsed into a single row
— keeping the most significant one (smallest adjusted p-value, ties broken
by the largest fold-change) and discarding the rest. The same collapsing
applies to any duplicate gene identifiers left after mapping, so a gene is
only ever represented once — a duplicate would otherwise be double-counted
by any downstream enrichment test.

RESULTS
-----------
The resolved table (downloadable as CSV) has one row per gene: input
identifier, gene symbol, Entrez ID, Ensembl gene ID, log2 fold-change,
p-value, adjusted p-value, and mapping status.

Summary cards report how many genes were mapped, unmapped, collapsed
(duplicate genes and transcripts), and how many rows were dropped for a
non-numeric fold-change or p-value. If no adjusted p-value column was
found in your upload, a note confirms one was computed for you.

This step produces no charts — it is a data-preparation stage, and the
summary cards ARE the result you should read. Check them before moving on,
because every downstream analysis inherits whatever happens here:

  - A high UNMAPPED count means those genes are absent from every analysis
    that follows. They were not silently discarded — they are still in the
    table with an UNMAPPED status so you can see exactly which ones — but
    nothing downstream can test a gene it could not identify.
  - A high COLLAPSED count is normal for transcript-level input and is
    what you want: without it, one gene represented by five transcripts
    would be counted five times by every enrichment test.
  - Rows DROPPED for a non-numeric fold-change or p-value are worth
    investigating if the number is more than a handful — it usually means a
    column was mis-detected, or the file uses a decimal comma or an "NA"
    string the parser read as text.
  - If the card confirms an adjusted p-value was COMPUTED, your upload
    didn't contain one. That is fine, but be aware the padj values in the
    resolved table are Omiscera's, not your original pipeline's, and were
    computed across the genes that survived this step.

THINGS TO KEEP IN MIND
---------------------------
  - A row is dropped only when its fold-change or p-value can't be read as
    a number, or its identifier cell is blank — every other row survives,
    including unmapped ones.
  - Less-studied, non-model organisms will generally have sparser
    identifier coverage, so expect a higher UNMAPPED count for those
    species — this reflects real, uneven database coverage, not a flaw in
    the mapping.
