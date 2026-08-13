GENOMICS — GENE FUNCTION ANNOTATION
=======================================

WHAT THIS MODULE IS FOR
--------------------------
This module asks a simple, per-gene question:

    "What does THIS gene actually do?"

For each gene, it looks up its Gene Ontology (GO) terms — Biological
Process, Molecular Function, Cellular Component — plus its UniProt
functional description, protein name, and keywords. There is no p-value,
no enrichment ratio, and no background universe: this is a direct lookup
of what's already known about each gene, not a statistical test. Useful
for inspecting individual genes of interest, or simply looking a gene up.

USER INPUT
-------------
Two ways to provide genes:

  - TABLE MODE: a table with one required column — Gene Symbol or Ensembl
    Gene ID — and one optional column, UniProt ID. Select an organism
    (needed to disambiguate a Gene Symbol across species), upload your
    table, and each gene is matched automatically, with any ambiguous
    matches flagged for you to resolve.

  - SINGLE-GENE MODE: type one identifier directly instead of uploading a
    table. The tool detects what you typed:
      - A UniProt accession (e.g. P38398) — no organism needed.
      - An Ensembl Gene ID (e.g. ENSG00000012048) — no organism needed,
        since these IDs are species-specific already.
      - Anything else is treated as a Gene Symbol, resolved via the
        organism you select (required, since the same symbol can exist in
        many species).

THE "DARK GENE" INDICATOR — ANNOTATION STATUS
---------------------------------------------------
GO annotation coverage is famously uneven: a small fraction of genes
account for most of the accumulated evidence, because well-studied genes
keep gathering more experimental support while many genes remain barely
characterized. A gene list you upload is therefore usually a MIX of
thoroughly studied genes and "dark," under-studied ones — a flat table of
GO terms alone doesn't tell you which is which.

Every gene is classified into one of three tiers:
  - WELL-CHARACTERIZED — has at least one experimental GO annotation (real
    lab evidence, not just computational inference) AND a curated function
    description. The "solid ground" tier.
  - MODERATE — some non-electronic support exists (a curated-but-not-
    fresh-experiment annotation, or a function description) but not enough
    to reach the well-characterized bar.
  - SPARSE (DARK) — only computationally-inferred annotations and/or no
    function description. These are the genes current databases know the
    LEAST about — worth flagging precisely because they're easy to
    overlook in a long results table, and are the most interesting
    candidates for further investigation.

The results table color-codes each row by this status. The Annotation
Status and GO-term counts are always computed from the gene's COMPLETE
evidence, regardless of the "Include IEA" display setting below — a gene's
characterization tier doesn't change just because you're choosing to hide
some term names from the table.

THE "INCLUDE IEA" SETTING
------------------------------
IEA ("Inferred from Electronic Annotation") marks annotations assigned
automatically by computational pipelines, without a human curator's
review — most GO annotations are IEA. They give broad coverage but are
weaker evidence than a curator-reviewed or experimental annotation.

Every gene's full evidence is always looked up. "Include IEA" is a
checkbox that controls only what's DISPLAYED in the GO term columns —
turning it off instantly re-filters the table to non-IEA evidence only.
The downloadable CSV always contains all evidence, at the term level,
regardless of this setting.

RESULTS AND DOWNLOADS
-------------------------
  - Summary cards — how many genes were annotated, and how many fall into
    each Annotation Status tier. Read these as a shape check on your gene
    list before anything else: a list that is mostly "Sparse (dark)" will
    produce thin results in every enrichment analysis you run downstream,
    because there is little annotation for those analyses to find. That is
    a property of what is currently known about your genes, not a fault of
    your experiment.
  - Gene cards (one per gene) — Gene, UniProt ID, Protein Name, Organism,
    GO Biological Process / Molecular Function / Cellular Component terms
    (filtered by the Include IEA setting), Keywords, the curated UniProt
    Function description, Annotation Status, and GO-term counts. Cards are
    colour-coded by Annotation Status so the dark genes are visible at a
    glance rather than buried in a long table.
  - Annotation table (CSV) — the same per-gene information in one row per
    gene, for spreadsheet use.
  - Tidy term-level table (CSV) — one row per gene-GO-term pair, covering
    EVERY annotation regardless of the Include IEA setting. Use this if you
    want to filter by evidence code, re-aggregate, or cite specific
    evidence yourself. This is the download to take if you intend to do
    your own analysis on the annotations rather than just read them.

No charts are produced here. This module answers a per-gene question ("what
is known about this gene?"), and a per-gene lookup has no group-level
pattern to plot — the colour-coded cards and the status counts are the
visual summary.

THINGS TO KEEP IN MIND
---------------------------
  - "Sparse (dark)" doesn't mean a gene is unimportant — only that current
    public databases have little direct evidence for it yet.
  - A negated GO annotation (GO's way of recording that a gene does NOT
    have a function some evidence suggested) is excluded from every term
    list and from the Annotation Status counts — including it as positive
    would be wrong, not just noisy.
  - Less-studied, non-model organisms will generally have sparser coverage,
    so more of their genes will land in the Moderate or Sparse tiers — this
    reflects real, uneven database coverage, not a flaw in the
    classification.
  - Single-gene UniProt-ID mode is exact, with no ambiguity possible.
    Single-gene Symbol mode can occasionally be ambiguous (a symbol
    matching multiple entries in that organism) — if so, look up the
    correct UniProt accession directly instead.
