TRANSCRIPTOMICS — OVER-REPRESENTATION ANALYSIS (ORA)
========================================================

WHAT THIS MODULE IS FOR
--------------------------
You give Omiscera a differential-expression table — every gene you measured,
with its fold-change and p-value. ORA first applies your padj/log2FC
thresholds to decide which genes count as "differentially expressed" (DEGs),
then asks:

    "Is there a common biological theme among my DEGs — are they made up of
    genes involved in a specific process, pathway, or gene set more than I'd
    expect if I'd picked genes at random from everything I actually tested?"

This is repeated for every term in each selected database, separately for
all DEGs, up-regulated DEGs, and down-regulated DEGs, producing a ranked
table of over-represented biological themes per direction.

Unlike GSEA (see the GSEA README), ORA only ever looks at genes that passed
your significance/fold-change cutoff — a coordinated but individually
sub-threshold shift across a whole pathway is invisible to it. Run both if
you're unsure which pattern your data has.

It runs against three curated knowledge bases:
  - Gene Ontology (GO) — standardized descriptions of what a gene product
    does, where it acts, and what process it participates in.
  - Reactome — curated metabolic and signaling pathways.
  - MSigDB — 8 curated collections you pick from (Hallmark, Perturbations,
    WikiPathways, PID, Phenotype, Oncogenic, Immunologic, Cell type),
    covering human and (mostly) mouse. See "MSIGDB COLLECTIONS AND
    LICENSING" below for exactly what each one is, what was deliberately
    left out, and why.

You can also upload your own custom gene-set database instead of/alongside
these three (see "CUSTOM DATABASE" below).

THE UNIVERSE (BACKGROUND)
------------------------------
The universe is every gene actually tested in this experiment (every mapped
row of your resolved table), NOT a generic genome-wide background. A gene
never expressed/detected in your experiment could not have been called DE
regardless of its biology, so it doesn't belong in the background either -
this is the standard RNA-seq convention (goseq, clusterProfiler both
recommend it), and it materially changes results versus a genome-wide
background. An optional uploaded background list (see below) narrows this
further, e.g. to a targeted gene panel.

THE STATISTICAL TEST
-----------------------
A one-sided Fisher's exact test (hypergeometric), same core test as
Genomics ORA. RNA-seq has a well-documented extra wrinkle GWAS-style gene
lists don't: longer genes accumulate more reads and are more likely to be
called DE independent of their biology, which biases a plain hypergeometric
test toward gene sets enriched for long genes. "Correct for gene-length
bias" (on by default) switches to a Wallenius non-central hypergeometric
test instead, GOseq-style: gene lengths are fetched (Ensembl BioMart) and a
length-vs-detectability trend is fit per direction (all/up/down each get
their own fit, since they're different gene lists with different detection
power). When no meaningful length bias is detected for a given direction, or
length data isn't available for the organism, that direction automatically
falls back to the plain hypergeometric test - reported after the run, not a
silent substitution.

FOLD ENRICHMENT
------------------
Alongside the p-value, each term gets a fold enrichment: how many times more
often genes from this term appear among your DEGs than random chance (given
the universe) would predict. Fold enrichment measures how STRONG the
enrichment is; the p-value measures how CONFIDENT you can be it isn't a
coincidence - complementary, not interchangeable.

MSIGDB COLLECTIONS AND LICENSING
-------------------------------------
MSigDB's main collection is Creative Commons Attribution 4.0 (CC BY 4.0),
copyright (c) Broad Institute, Inc., Massachusetts Institute of Technology,
and the Regents of the University of California — free for academic and
commercial use, redistributable, attribution required. This app fetches it
directly from MSigDB's own release server (verified openly accessible, not
registration-gated), version PINNED at 2026.1 so results stay reproducible
run to run — the pin is bumped deliberately, not automatically.

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

Has no effect on evidence-code filtering (see below): MSigDB carries no
IEA/evidence-code concept at all, unlike GO and Reactome.

EVIDENCE CODES (IEA)
-------------------------
GO and Reactome both distinguish curator-reviewed annotations from
"IEA" (Inferred from Electronic Annotation) ones assigned automatically,
without review. Included by default; uncheck to restrict both to curated
evidence only. Reactome's curated (non-IEA) coverage is essentially
human-only, so unchecking this for another organism will likely leave
Reactome empty. No effect on MSigDB or a custom database.

BACKGROUND (REFERENCE) GENE LIST
-------------------------------------
Optional. The default universe (genes tested in this experiment) is already
the RNA-seq-appropriate background, so this is only needed to narrow
further - e.g. to a targeted panel you're specifically interested in.
Applied globally, once, across every selected database (unlike Genomics,
where background is set per-database) - this module already has one shared
universe, so there is one shared narrowing on top of it.

CUSTOM DATABASE
--------------------
Upload your own gene sets (GMT format, or a tidy term/gene table), always in
gene-symbol space. The background universe becomes the union of every gene
in the file, intersected with the genes tested in this experiment - same
convention as GO/Reactome/MSigDB.

CORRECTING FOR TESTING MANY TERMS AT ONCE
-----------------------------------------------
Benjamini-Hochberg FDR is applied separately within each independent
family - every database, crossed with direction (all/up/down), and (since
MSigDB was added) every SELECTED MSIGDB COLLECTION, crossed with direction,
is its own family with its own correction. Selecting several MSigDB
collections multiplies how many such families one run can produce: up to 24
(8 collections x 3 directions), on top of GO/Reactome/Custom database. The
results screen reports exactly how many independently-corrected tables a
given run produced. Scanning every table and reporting the single best hit
across all of them inflates the true error rate above what any one FDR
column states - worth remembering before selecting every collection "just
to see." A Bonferroni column is also reported, stricter but more
conservative than FDR.

RESULTS AND VISUALIZATION
------------------------------
Each database/direction produces a results table (downloadable as CSV).
Five figure types are available, shared with Genomics ORA and GSEA so
results look and export consistently across the app:
  - Dot plot / bar plot - one row per term, encoding fold enrichment, FDR,
    and gene count. The usual starting point.
  - Volcano plot - fold enrichment vs. significance (FDR), for a quick view
    of which terms are both strongly and confidently enriched.
  - Enrichment map - links terms that share DE genes, revealing clusters of
    related, overlapping biology rather than independent findings.
  - Gene-concept network - shows which specific genes are driving each
    significant term.

MSigDB collection summary - shown once per run (not per database), only
when MSigDB was selected: one horizontal bar per selected collection, its
length showing how many gene sets cleared your significance threshold,
colour showing the direction of its best hit (red = positive fold
enrichment / up-regulated, blue = negative). Exists because with up to 8
collections in one run, nothing else answers "which collection is actually
carrying the signal" without clicking through every one by hand.

Circle SIZE on the dot plot/enrichment map/gene network encodes gene count;
circle AREA (not diameter) tracks it, so a circle with four times the area
represents four times the genes. Hover any point for the exact number.

THINGS TO KEEP IN MIND
---------------------------
  - ORA only sees genes that passed your padj/log2FC cutoff - a small,
    coordinated shift across a whole pathway with no individual gene
    crossing the threshold is invisible to it. Try GSEA if that's a concern.
  - MSigDB covers human and mouse only, and not every collection has a
    mouse file - see "MSIGDB COLLECTIONS AND LICENSING" above for the exact
    per-collection availability table and the HPO-vs-MPT phenotype caveat.
  - GO's gene-set counts are inflated by ontology propagation (a gene
    annotated to a specific term also counts toward every broader parent
    term), so don't be surprised to see closely related parent/child GO
    terms both come out significant with overlapping gene lists.
