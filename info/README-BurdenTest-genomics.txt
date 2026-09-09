GENOMICS — GENE-SET BURDEN TEST
==================================

WHAT THIS MODULE IS FOR
--------------------------
Sometimes what you have is not a clean list of genes that are "in" or "out"
— it's a CONTINUOUS score for every gene: a variant severity score (e.g.
CADD), a GWAS association strength (e.g. a p-value), or any other weighted
measure of how strongly a gene is implicated.

The Gene-set Burden Test answers this question:

    "Do the genes in this biological pathway carry systematically HIGHER
    scores than the rest of my scored genes — more variant burden, stronger
    association — than expected if the scores had nothing to do with
    pathway membership?"

Nothing is thresholded first: every scored gene contributes, and a pathway
can come out significant because all of its genes are moderately elevated,
not only because a few are extreme. This "competitive" framing — comparing a
pathway's score distribution against the rest of the scored genes — is the
same one used by established gene-set tools such as CAMERA and MAGMA.

It runs against the same curated knowledge bases available elsewhere in
Genomics: Gene Ontology, Reactome and KEGG PATHWAY, or your own uploaded
gene-set database. (GWAS Catalog is not offered here — it is symbol-native
and an unconventional unit for a burden test.)

USER INPUT
-------------
A table with:
  - One identifier column (required): Gene Symbol OR Ensembl Gene ID.
  - One numeric column (required): your Variant Score — any per-gene metric
    (a CADD score, a GWAS p-value, a weighted variant count, or similar).
  - UniProt ID (optional) — if you already know it, supplying it skips
    resolution for that row.

Select an organism, upload your table, and each gene is matched
automatically to a UniProt entry, with any ambiguous matches flagged for you
to confirm. There is no Homology (BLAST) option here, since a raw sequence
isn't a stable identifier to attach a score to.

CHOOSING THE SCORE DIRECTION — "HIGHER" OR "LOWER" = MORE BURDEN
---------------------------------------------------------------------
The tool reads your score column as a plain number — it does not try to
guess what kind of metric it is. What tells it how to interpret the number
is the DIRECTION setting you choose next to the column:

  - Pick "Higher = more burden" for scores where a BIGGER number means more
    signal: a CADD score, a weighted variant count, a -log10(p-value).
  - Pick "Lower = more burden" for scores where a SMALLER number means more
    signal: a raw, untransformed p-value (e.g. straight from a GWAS).

You do not need to convert a raw p-value to -log10(p) yourself — choosing
"Lower = more burden" on the raw p-value gives the same significant terms
as pre-transforming it and choosing "Higher." (The Effect column's exact
number will differ between the two, but which terms come out significant
does not.)

If a result looks "backwards" — a pathway you expected to show high burden
instead shows a negative effect, or vice versa — check this setting first.
It is the single most common source of a flipped result.

Rows where the score can't be read as a number (blank, text, etc.) are
skipped and reported to you as "non-numeric" rows.

THE STATISTICAL TEST
-----------------------
For every pathway/term, the test compares the scores of genes INSIDE it
against the scores of every OTHER scored gene, using a one-sided Wilcoxon
rank-sum test (Mann-Whitney U). This test:
  - Makes no assumption about how your scores are distributed.
  - Always gives the same answer for the same input (no randomness
    involved).
  - Is the standard choice for this kind of comparison, also used by
    GSEA/CAMERA/MAGMA-style tools.

Each term's result reports how many of your scored genes belong to it, the
average score inside vs. outside the term (the difference is the
"Effect"), a p-value, and an FDR-corrected p-value.

MULTIPLE-TESTING CORRECTION
-------------------------------
Because many pathways are tested at once, some will look "significant"
purely by chance. Benjamini-Hochberg FDR correction is applied separately
within each database (GO:BP, GO:MF, GO:CC, Reactome, KEGG, or your custom
database each get their own correction). A Bonferroni option is also
available for a stricter, more conservative threshold. Significance
defaults to FDR < 0.05, adjustable per database.

KEGG identifies organisms by its own organism code (e.g. "hsa"), a separate
choice from the UniProt organism above, made during retrieval rather than
here — see "USING KEGG" below. KEGG's coverage of non-model organisms is
usually the best of the three databases here, the same reason it is offered
in Metabolic Connection and Genomics ORA.

USING KEGG
--------------
KEGG tests in its own gene-ID space, not UniProt accession space (an earlier
version matched via UniProt and inflated term sizes wherever a gene had more
than one UniProt entry, and could never reach genes with no UniProt
cross-reference — worst for non-model organisms). To use it: resolve your
genes normally first, then use the "Resolve KEGG gene IDs" box that appears
on the retrieval screen — pick a KEGG organism and click "Resolve KEGG gene
IDs." This is optional; skip it if you are not using KEGG here.

This step can also rescue genes that came back "ambiguous" under UniProt
alone — common for non-model organisms, which usually have no UniProt
reference proteome to break a tie between candidates. If KEGG's own data
identifies the gene uniquely, it is resolved (shown as "Resolved via KEGG")
and its score becomes usable everywhere, not only for KEGG. KEGG never
guesses: a gene stays ambiguous if more than one candidate remains plausible
even with KEGG's help. Results still show your own gene symbol or Ensembl ID
in the "Genes" column, never a raw KEGG gene ID.

TERM-SIZE FILTERING
-----------------------
A minimum term size (default 3 genes) excludes tiny, unreliable sets; an
optional maximum lets you exclude very large, generic sets. A term also
needs at least one gene outside it to serve as the comparison group.

RESULTS AND VISUALIZATION
-----------------------------
Each database produces a results table (downloadable as CSV) plus five
figures, grouped on screen the same way they are grouped below.

1. GENE SET OVERVIEW — "which pathways carry the burden, and how strongly?"
  - Dot plot — one row per term. Horizontal position is the Effect (the
    difference between the term's mean score and everything else's), colour
    is significance, and circle size is how many scored genes the term
    contains. A term far to the right with a dark, large circle is one where
    many genes each carry a clearly elevated score.
  - Bar plot — the same terms as plain bars, rankable by FDR or by Effect.
    Ranking by Effect answers "where is the burden largest?"; ranking by FDR
    answers "where am I most confident?" — these are different questions and
    often give a different top of the list.
  - Volcano plot — Effect against significance for every term at once. The
    top-right corner holds terms that are both strong and confident. A term
    high but far left is confidently a *small* shift; one far right but low
    is a large shift resting on few genes.

2. GENE SET RELATIONSHIPS — "are these separate findings, or one story?"
  - Enrichment map — draws a line between two high-burden terms when they
    share genes. A tight cluster is usually ONE underlying signal reported
    under several names, not several independent results.
  - Gene-concept network — shows the individual genes behind each term,
    with each gene node coloured by its OWN burden score. This is where you
    see whether a term is significant because many genes are moderately
    elevated, or because one or two extreme genes are dragging the whole
    set — a distinction the table alone cannot show you, and one that
    changes how much you should trust the result.

All charts are interactive (draggable nodes, adjustable font size and number
of terms shown, adjustable significance threshold) and exportable as SVG,
PNG, or 300 dpi TIFF.

On the dot plot, enrichment map and gene-concept network, circle SIZE
encodes how many scored genes are in the term and circle COLOUR encodes
significance. Size is scaled so the largest circle in a given chart is a
fixed, readable size no matter how many genes you uploaded, and so that
circle AREA (not diameter) tracks the count. Hover any point for the exact
number.

A REMINDER ABOUT SIGN: if you set "lower = more burden" (e.g. a raw
p-value), the scores were negated internally so that higher always means
more burden. The Mean score columns in the table therefore show negated
values, not your original numbers — the ranking and the Effect sign are
still correct, but don't read those two columns as raw input values.

THINGS TO KEEP IN MIND
---------------------------
  - Effect (a plain average difference) and FDR (from a rank-based test)
    usually agree but aren't identical — read Effect as "how big" and FDR
    as "how confident."
  - The test treats each gene's score as independent; genes that are
    correlated with each other (e.g. co-regulated) can make results
    slightly more likely to look significant than they truly are — a
    normal, standard trade-off for this kind of fast test.
  - Works best with a large, representative set of scored genes — ideally
    every gene your assay could plausibly have measured, not a hand-picked
    shortlist.
  - A significant, positive-effect pathway means association, not proof of
    causality.
  - Double-check your score direction setting (see above) if a result looks
    backwards — it is by far the most common cause.
  - KEGG's REST API and academic use — see README-ORA-genomics.txt; the same
    points apply here unchanged.
