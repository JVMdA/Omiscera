GENOMICS — OVER-REPRESENTATION ANALYSIS (ORA)
================================================

WHAT THIS MODULE IS FOR
--------------------------
You give Omiscera a list of genes — genes that were differentially
expressed in an experiment, came out of a GWAS hit list, were targets of a
drug, or any other list you assembled. On its own, a list of genes doesn't
tell you much. The question ORA answers is:

    "Is there a common biological theme in this list — is it made up of
    genes involved in a specific process, pathway, or trait more than I'd
    expect if I'd picked genes at random?"

It answers this against three curated knowledge bases:
  - Gene Ontology (GO) — standardized descriptions of what a gene product
    does, where it acts, and what process it participates in.
  - Reactome — curated metabolic and signaling pathways.
  - GWAS Catalog — genes associated with human traits and diseases in
    published genome-wide association studies.

You can also upload your own custom gene-set database instead of these
three (see "CUSTOM DATABASE" below).

THE SCIENTIFIC QUESTION
----------------------------
For every biological term in a database (e.g. the GO term "oxidative
phosphorylation"), ORA asks: "Of the genes in my list, how many belong to
this term — and is that more than chance alone would predict, given the
term's size and my list's size?" This is repeated for every term in the
database, producing a ranked table of which biological themes are
over-represented in your list.

THE STATISTICAL TEST
-----------------------
ORA uses a one-sided Fisher's exact test. For each term, it compares how
many of your genes belong to it against how many would be expected by
chance, given the term's size and the size of the background ("universe,"
see below). A small p-value means the overlap is bigger than chance would
produce.

FOLD ENRICHMENT
------------------
Alongside the p-value, each term gets a fold enrichment: how many times
more often genes from this term appear in your list than random chance
would predict. A fold enrichment of 3 means three times more often than
expected. Fold enrichment measures how STRONG the enrichment is; the
p-value measures how CONFIDENT you can be it isn't a coincidence — they're
complementary, not interchangeable. A term can have a huge fold enrichment
but a weak p-value if only 1-2 genes are involved, or a modest fold
enrichment with a very strong p-value if many genes are involved.

CORRECTING FOR TESTING MANY TERMS AT ONCE
-----------------------------------------------
Because ORA tests every term in a database at once (GO alone has tens of
thousands), some will look "significant" by chance alone. Omiscera corrects
for this with:
  - Benjamini-Hochberg FDR (default) — controls the expected proportion of
    false positives among the terms you call significant, applied
    separately within each term family (GO:BP, GO:MF, GO:CC, Reactome, GWAS
    Catalog, or your custom database).
  - Bonferroni — a stricter, more conservative alternative, also reported.

By default, FDR < 0.05 is treated as significant; adjustable per database.

BACKGROUND / UNIVERSE
------------------------
The "universe" is the pool of genes treated as "everything that could have
been in your list." By default, this is the whole annotated genome for
your species.

This default is fine when your list could plausibly have contained any
gene (e.g. a curated disease-gene list from the literature). It is NOT
appropriate when your list came from a restricted pool to begin with — the
classic example is RNA-seq: if your experiment only detected 14,000 of
~20,000 genes, testing against the whole genome makes anything relevant to
your tissue look artificially "enriched," simply because those genes were
more detectable in your assay, not because your condition favored them.

You can upload your own background/reference gene list (for GO and
Reactome) to fix this — the universe then narrows to the genes your
experiment could actually have observed. This is optional and off by
default.

(Background upload applies to GO and Reactome only. GWAS Catalog works in
gene-symbol space natively, and a custom database always defines its own
universe from the file you upload.)

TERM-SIZE FILTERING
-----------------------
Very small terms (1-2 genes) give statistically unreliable results — one
coincidental match can look significant. Very large, generic terms (e.g.
"protein binding") are easy to "hit" but tell you little. By default,
terms smaller than 3 genes are excluded; you can also set an optional
maximum to filter out overly generic categories.

SPECIES SELECTION
--------------------
Choosing your organism tells GO, Reactome, and GWAS Catalog which species'
data to use, and is also used to match your uploaded genes to the right
species when the same gene symbol exists in multiple organisms.

RESULTS AND VISUALIZATION
-----------------------------
Each database produces a results table (Term, gene counts, fold enrichment,
p-value, FDR, overlapping genes — downloadable as CSV), plus five figures,
grouped on screen the same way they are grouped below.

1. GENE SET OVERVIEW — "which terms came out, and how strongly?"
  - Dot plot — one row per term. Colour shows significance (darker = smaller
    FDR) and circle size shows how many of your genes are in the term. The
    horizontal axis follows your "Rank & show by" choice — with one
    deliberate exception: when you rank by FDR, the axis shows gene ratio
    instead, because colour is already carrying FDR and plotting it twice
    would waste the axis. The axis is always labelled with what it is
    showing. The usual starting point: a large, dark circle far to the
    right is a term that is both strongly and confidently enriched.
  - Bar plot — the same terms as plain bars. Easier to read exact values
    off, and the better choice for a figure in a paper where circle size
    would need its own explanation.
  - Volcano plot — fold enrichment (how strong) against significance (how
    confident) for every term at once. The top-right corner holds terms
    that are both. This is the fastest way to spot the trap described under
    FOLD ENRICHMENT above: terms high up but far left are confident but
    weak, terms far right but low are strong but shaky.

2. GENE SET RELATIONSHIPS — "are these separate findings, or one story?"
  - Enrichment map — draws a line between two terms when they share genes.
    Terms that clump together are largely the same biology reported under
    several names; a tight cluster of ten terms is usually ONE finding, not
    ten. This is how you avoid over-counting your results. The "Edge
    threshold" control sets how many shared genes are needed before a line
    is drawn — raise it to keep only the strongest overlaps.
  - Gene-concept network — shows the individual genes behind each enriched
    term, and which genes are shared between terms. Use it to find the few
    genes doing most of the work across several terms; those are usually
    the ones worth following up experimentally.

All charts are interactive (draggable nodes, adjustable font size and number
of terms shown, adjustable significance threshold) and exportable as SVG,
PNG, or 300 dpi TIFF.

On the dot plot, enrichment map and gene-concept network, circle SIZE
encodes gene count and circle COLOUR encodes significance. Size is scaled so
the largest circle in a given chart is a fixed, readable size no matter how
large your gene list is, and so that circle AREA (not diameter) tracks the
count — a circle with four times the area represents four times the genes.
Hover any point for the exact number.

THINGS TO KEEP IN MIND
---------------------------
  - A low p-value isn't the same as a large effect — always check fold
    enrichment alongside FDR.
  - Your choice of background changes the answer. This isn't a flaw — it's
    what "over-represented relative to chance" means. Pick the background
    that matches how your gene list was actually generated.
  - FDR is corrected per term family, not globally.
  - Non-human, less-studied organisms will have sparser annotation,
    especially for GWAS Catalog (essentially human-only).
