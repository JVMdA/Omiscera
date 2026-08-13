TRANSCRIPTOMICS — GENE SET ENRICHMENT ANALYSIS (GSEA)
=========================================================

WHAT THIS MODULE IS FOR
--------------------------
You give Omiscera a differential-expression table — every gene you measured,
with its fold-change and p-value. GSEA uses ALL of them, with no significance
cutoff applied first. The question it answers is:

    "Do the genes in this pathway tend to cluster near the TOP or BOTTOM of
    my whole ranked gene list — even if none of them individually reached
    statistical significance?"

Every resolved gene is ranked from most up- to most down-regulated (or by
whichever metric you choose — see below), and each pathway is scored by
where its members fall in that ranking. Nothing is discarded before testing.

This makes GSEA well suited to the common case where a whole pathway shifts
by a small, coordinated amount — every gene moving a little in the same
direction — which is a real biological signal that a per-gene significance
cutoff would discard gene by gene, since no single member stands out on its
own.

It runs against three curated knowledge bases:
  - Gene Ontology (GO) — standardized descriptions of what a gene product
    does, where it acts, and what process it participates in.
  - Reactome — curated metabolic and signaling pathways.
  - MSigDB — 8 curated collections you pick from (Hallmark, Perturbations,
    WikiPathways, PID, Phenotype, Oncogenic, Immunologic, Cell type),
    covering human and (mostly) mouse. See "MSIGDB COLLECTIONS AND
    LICENSING" below for exactly what each one is, what was deliberately
    left out, and why.

You can also upload your own custom gene-set database.

MSIGDB COLLECTIONS AND LICENSING
-------------------------------------
MSigDB's main collection is Creative Commons Attribution 4.0 (CC BY 4.0),
copyright (c) Broad Institute, Inc., Massachusetts Institute of Technology,
and the Regents of the University of California — free for academic and
commercial use, redistributable, attribution required. This app fetches it
directly from MSigDB's own release server (verified openly accessible, not
registration-gated), version PINNED at 2026.1 so results stay reproducible
run to run — the pin is bumped deliberately, not automatically, and every
result footer/CSV records which version produced it.

Offered (8 of MSigDB's collections). Counts are for release 2026.1 and were
measured directly from the downloaded files, not copied from a webpage
(tests/live_msigdb_urls.py re-checks them on every run):

  Collection              Human file         Mouse file          Sets (Hs/Mm)
  Hallmark                h.all              mh.all                 50 / 50
  Perturbations (CGP)     c2.cgp             m2.cgp               3555 / 984
  WikiPathways            c2.cp.wikipathways m2.cp.wikipathways    925 / 202
  PID                     c2.cp.pid          (none)                 196 / —
  Phenotype (HPO / MPT)   c5.hpo             m5.mpt                5793 / 92
  Oncogenic               c6.all             (none)                 189 / —
  Immunologic             c7.all             m7.all               5219 / 787
  Cell type               c8.all             m8.all                 866 / 233

Only PID and Oncogenic have no mouse file at all ("—"). Every other
collection runs for both organisms - though note how much smaller several of
the mouse sets are, which is a real coverage difference, not a fetch problem.

A note on "Phenotype": the human file is HPO (Human Phenotype Ontology) and
the mouse file is MPT (Mammalian Phenotype Ontology) — the same KIND of
information (clinical/organismal phenotype terms), genuinely distinct from
GO, but two DIFFERENT ontologies. A human run and a mouse run of "Phenotype"
are not testing the same term set, only an analogous one.

Mouse coverage is NOT a mirror of human's — PID and Oncogenic have no mouse
file, and several others are much smaller for mouse, verified against
MSigDB's own collections pages rather than assumed from the file naming
(mouse collection numbers do not track human ones: mouse's cell-type
collection is M8, matching human C8, but mouse M3 is regulatory targets, not
an analog of anything offered here). The app greys out exactly the
collections with no file for your organism; it never guesses a filename.

Deliberately excluded, and why:
  - GO (10,490 sets) and Reactome's own MSigDB copy (1,839 sets) duplicate
    this app's own GO/Reactome integrations, from a different source and
    version — they would give DIFFERENT answers sitting right next to the
    "real" GO/Reactome families, which would read as a bug, not a feature.
  - KEGG_LEGACY, KEGG_MEDICUS, and BioCarta are excluded on licensing
    grounds specifically: KEGG_LEGACY and BioCarta are © their original
    authors under "qualified permission" only, not CC BY; KEGG_MEDICUS is
    CC BY-SA (a ShareAlike/viral term, different from everything else
    offered here). KEGG is also excluded from this app on principle.
  - TF-target sets (C3:TFT, 1,116 sets) are left out because this app
    already has a dedicated Transcription Factor / Regulon module with its
    own directional statistics and empirical-null calibration — a second,
    unsigned TF-target source sitting in ORA/GSEA would produce a competing
    TF answer with no reconciliation between the two, a real interpretation
    hazard rather than a convenience.
  - Positional (C1) and Computational (C4) are niche (CNV/karyotype-driven
    studies, and cancer-specific/partly-legacy sets respectively).

MULTIPLE-TESTING FAMILIES, SHARPENED BY MSIGDB
---------------------------------------------------
Every database/collection here gets its OWN independent Benjamini-Hochberg
correction — there is no correction ACROSS them (same convention as ORA and
the Regulon module). Selecting several MSigDB collections multiplies how
many such independent families one run can produce: up to 8 (one per
collection), on top of GO/Reactome/Custom database. The results screen
reports exactly how many independently-corrected tables a given run
produced. Scanning every table and reporting the single best hit across all
of them inflates the true error rate above what any one FDR column states —
worth remembering before selecting every collection "just to see."

AN IMPORTANT LIMITATION, STATED UP FRONT
--------------------------------------------
The classic version of GSEA was designed for experiments with individual
SAMPLES per condition (e.g. 5 treated vs. 5 control), and gets its
statistical null distribution by reshuffling which samples belong to which
group ("phenotype permutation"). This app's input is a differential-
expression SUMMARY table (one row per gene, with a fold-change and a
p-value) — there are no individual samples to reshuffle.

This module is therefore "pre-ranked" GSEA, the same approach used by
widely-used tools such as fgsea and clusterProfiler's GSEA function when
you only have a ranked gene list. Its null distribution comes from randomly
reshuffling which genes fall in a pathway ("gene-set permutation") instead.

This is a real, accepted limitation, not a shortcut specific to this app:
gene-set permutation assumes genes behave independently, so a pathway whose
genes are strongly co-regulated (move together for reasons unrelated to
your experiment) can look more significant than it really is. Treat a
borderline result on a tightly co-regulated pathway with a bit of extra
skepticism.

CHOOSING A RANKING METRIC
------------------------------
Every gene needs a single number to be ranked by. Four options are offered:

  - Signed -log10(p-value) (recommended default) — combines direction
    (from log2 fold-change) with statistical confidence (from the
    p-value). A gene that changed a lot AND was measured confidently ends
    up near the extreme of the ranking; a gene with a large fold-change but
    weak/noisy evidence does not dominate the ranking the way it would
    under a fold-change-only metric.
  - log2 Fold-Change — the raw effect size, ignoring how confidently it was
    measured. Simple, but a single noisy gene with a huge fold-change and a
    poor p-value can sit at the very top of the list.
  - π-value (log2FC × -log10 p) — rewards genes that are both strongly
    changed and confidently measured. Useful, but its scale is less
    intuitive and it can exaggerate the gap between your top genes and
    everything else.
  - Signed -log10(padj) — like the recommended default, but using the
    FDR-adjusted p-value. More conservative, but many genes commonly share
    an identical (or maximal) adjusted p-value, which flattens large
    stretches of the ranking and can weaken the test.

Genes with no usable value for the chosen metric (only relevant for the
padj-based option, since not every gene necessarily has an adjusted
p-value) are dropped from the ranking and the count is reported.

THE STATISTICS, IN PLAIN TERMS
-----------------------------------
For each pathway, GSEA walks down your ranked gene list from top to bottom,
keeping a running score that goes UP a little every time it passes one of
the pathway's genes, and DOWN a little every time it passes any other gene.
The pathway's Enrichment Score (ES) is the single largest value that
running score reaches in either direction — a large positive ES means the
pathway's genes are clustered near the top (up-regulated) end of your
ranking; a large negative ES means they're clustered near the bottom
(down-regulated) end.

Because bigger pathways naturally tend to reach bigger running-score
extremes, ES is normalized into NES (Normalized Enrichment Score) so that
pathways of different sizes can be compared fairly and ranked against each
other. NES's sign IS the direction: positive means the pathway is enriched
toward the up-regulated end of your ranking, negative toward the
down-regulated end. Both directions come out of a single run — there is no
"up"/"down" setting to choose.

Significance is assessed by comparison against a null: what ES values would
this same pathway size produce if pathway membership were assigned to genes
at random, many times over? This gives a nominal p-value and, pooled across
every tested pathway, an FDR q-value — the number you should generally rely
on for deciding what is "interesting."

THE LEADING EDGE
--------------------
For each pathway, the "leading edge" is the subset of its genes that
actually drove the enrichment score to its peak — the genes at or before
the running score's maximum (for an up-regulated pathway) or at or after
its minimum (for a down-regulated one). This is usually a much shorter,
more actionable list than the pathway's full gene membership, and is what
the enrichment map and gene network charts are built from.

MULTIPLE-TESTING CORRECTION
-------------------------------
Two numbers are reported per pathway:
  - FDR q-value — the primary number to use. It pools the permutation null
    across every tested pathway (not just the one being scored), which is
    what makes it comparable to results from other GSEA tools. The
    established convention for this method treats FDR q < 0.25 as worth a
    look. That threshold is deliberately looser than the 0.05 used for many
    other tests, because pooling the null across all pathways already makes
    this a fairly conservative correction on its own.
  - FDR (BH) — a plain Benjamini-Hochberg correction on the nominal
    p-value, shown alongside for comparison. With a limited number of
    permutations, many pathways' nominal p-values hit the same floor value,
    so this column is often less discriminating than the FDR q-value.
A Bonferroni column is deliberately not offered here: at the resolution a
practical number of permutations allows, it would be a column of
essentially meaningless 1.0's.

GENE-SET SIZE FILTERING
----------------------------
A minimum gene-set size (default 15) excludes very small pathways, whose
enrichment scores are statistically unstable — with only a handful of genes,
the running score can reach an extreme by coincidence. An optional maximum
excludes very large, generic ones. Only pathway genes actually present in
your ranked list count toward this size, so a 200-gene pathway of which your
experiment measured 12 counts as size 12 and is excluded at the default.

PERMUTATIONS, WEIGHT, AND THE RANDOM SEED
----------------------------------------------
  - Number of permutations controls how finely a p-value can be resolved —
    more permutations give finer resolution at the cost of a longer run.
  - Weight exponent (default 1.0) is the standard GSEA weighting, where
    each gene's contribution to the running score is proportional to the
    magnitude of its ranking metric. Setting it to 0 gives an unweighted
    test, equivalent to a Kolmogorov-Smirnov-style statistic.
  - The random seed makes a re-run with identical settings reproduce
    identical results exactly — useful for verifying a result, or for
    reporting exactly how it was produced.

RESULTS AND VISUALIZATION
------------------------------
Each database produces a results table (downloadable as CSV) plus eleven
figures, grouped on screen the same way they are grouped below. Read this
section next to the app: the screen carries only short labels saying which
control drives which figure, and the explanation of what each figure MEANS
lives here.

MSigDB collection summary — shown once per run (not per database), only
when MSigDB was selected, above the individual result tables: one horizontal
bar per selected collection, its length showing how many gene sets cleared
your significance threshold, colour showing the direction of its best hit
(red = positive NES/fold enrichment, blue = negative). Exists because with
up to 8 collections in one run, nothing else answers "which collection is
actually carrying the signal" without clicking through every one by hand.

1. GENE SET OVERVIEW — "which pathways came out, and how strongly?"
  - NES dot plot / bar plot — one row per pathway, encoding NES, FDR
    q-value, and leading-edge size. The usual starting point.
  - NES volcano plot — enrichment strength (NES) vs. significance (FDR
    q-value), for a quick view of which pathways are both strongly and
    confidently enriched.

2. GENE SET RELATIONSHIPS — "are these findings independent, or one story?"
  - Enrichment map — links pathways that share leading-edge genes,
    revealing clusters of related, overlapping biology. Highly overlapping
    pathways are not independent findings; this is how you see that.
  - Leading-edge gene network — shows which specific genes are driving each
    significant pathway.
  The "Edge threshold" control belongs to this group only: it sets how many
  shared genes two pathways need before a line is drawn between them.

3. PER GENE SET — one pathway at a time, chosen in the "Gene set" menu
  - The enrichment plot — the classic three-panel GSEA figure: the running
    enrichment score, where that pathway's genes fall in your ranked list
    (the "Hits" rug), and the ranking metric itself across the whole list.
    The peak of the running score is the pathway's ES; the dotted vertical
    line marks where it occurs.
  - Null distribution — the permutation null for pathways of this size,
    with your observed NES drawn on it. This is the figure that shows you,
    directly, how far out in the null your result actually sits, rather
    than asking you to trust a q-value on its own.
  Both figures follow the "Gene set" menu ONLY. Only the most significant
  pathways in each database have a stored curve (storing this for every
  tested pathway would be far too much data when thousands are tested).

4. COMPARE GENE SETS
  - Multi-set enrichment plot — 2 to 5 pathways' running-score curves on
    one figure, each with its own row in the shared Hits rug, for comparing
    where different pathways concentrate in the same ranking. Capped at 5
    because beyond that the curves stop being separable by eye. Its
    checkboxes drive this figure only — they do not change group 3.

5. LEADING EDGE — "which genes are actually doing the work?"
  - Ridgeline plot — one density curve per pathway, showing where that
    pathway's genes sit along the ranking metric. Red means positive NES
    (the pathway is enriched toward the top of the ranking), blue means
    negative NES; a colour key is drawn on the figure. A ridge shifted
    right sits in the up-regulated end, one shifted left in the
    down-regulated end. Its "Gene population" control switches between the
    leading-edge (core) genes and all set members.
  - Leading-edge heatmap — pathways as rows, genes as columns, coloured by
    the ranking metric. Blank cells mean the gene is not in that pathway —
    they are NOT a metric of zero. Its three controls ("Genes shown",
    "Gene selection", "Row order") change this figure only; they do not
    affect the ridgeline, which is why each control is tagged on screen.

6. DIAGNOSTICS — quality control, not a biological finding
  - NES vs. set size — whether NES normalization actually removed the
    dependence of enrichment score on how many genes a pathway has. The
    median line should be roughly flat across the size range. A strong
    trend here means results at one end of the size range should be read
    more cautiously.

All charts are interactive — draggable nodes, adjustable font size and
number of pathways shown, adjustable significance threshold — and exportable
as SVG, PNG, or 300 dpi TIFF.

On the dot plot, enrichment map and gene network, circle SIZE encodes how
many genes are involved and circle COLOUR encodes significance. Size is
scaled so the largest circle in a given chart is a fixed, readable size no
matter how big your gene list is, and so that circle AREA (not diameter)
tracks the gene count — a circle with four times the area represents four
times the genes. Hover any point for the exact number.

THINGS TO KEEP IN MIND
---------------------------
  - This is pre-ranked, gene-set-permutation GSEA (see the limitation
    above), because the input is a summary table with no per-sample
    expression data to reshuffle.
  - Your choice of ranking metric materially affects which pathways come
    out on top; if a result surprises you, try a second metric before
    concluding it's wrong.
  - A significant, positive-NES pathway means association with the
    up-regulated end of your ranking, not proof of causality.
  - MSigDB covers human and mouse only, and not every collection has a
    mouse file (see "MSIGDB COLLECTIONS AND LICENSING" above for the exact
    per-collection availability table and the HPO-vs-MPT phenotype caveat).
  - GO's gene-set counts are inflated by ontology propagation (a gene
    annotated to a specific term also counts toward every broader parent
    term), so don't be surprised to see closely related parent/child GO
    terms both come out significant with overlapping leading edges.
