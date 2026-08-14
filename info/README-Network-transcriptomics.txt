TRANSCRIPTOMICS — NETWORK / MODULE DETECTION
=====================================================

WHAT THIS MODULE IS FOR
--------------------------
Over-Representation Analysis, Gene Set Enrichment Analysis and TF/Regulon all
ask "which pre-defined gene set is over-represented in my results?" This
module asks two different questions, using the known protein-protein
interaction network (STRING) for your organism instead of a gene set:

    "Starting from my differentially expressed genes as 'seeds,' which
    proteins does their signal concentrate on as it spreads across the
    network — INCLUDING proteins I never uploaded?"

    "Do the significant proteins organize into densely-interconnected
    functional MODULES?"

The numerics behind this — the network parser, the signal-spreading step,
the statistical test, and the clustering — are the SAME engine the Genomics
"Network Propagation / Module Detection" module uses (see that module's own
README for a from-scratch walkthrough of the maths). Only one copy of that
engine exists, so the two tiles can never quietly drift into disagreeing
about how propagation itself works. Everything documented below is what is
DIFFERENT here, because a differential-expression table is not a curated
candidate list:

    Genomics                        Transcriptomics
    -------------------------------  ------------------------------------
    a short curated candidate list   every gene that was TESTED (thousands)
    UniProt accessions               gene symbol / Ensembl ID / Entrez ID
    every uploaded gene is a seed    the table must be THRESHOLDED first
    optional unsigned score          log2 fold-change (SIGNED) + p/padj

USER INPUT
-------------
Shares the same differential-expression retrieval step as Over-Representation
Analysis, Gene Set Enrichment Analysis and TF/Regulon (see the Retriever
README) — the resolved table is reused directly. No extra upload is needed.

Seeds are chosen from that table using the SAME padj/log2FC cutoffs and
"all/up/down" direction control as ORA and GSEA (Transcriptomics_module's own
select_degs function — the exact same call, so this module can never disagree
with ORA/GSEA about what counts as a differentially expressed gene). A
DIRECTION choice restricts WHICH genes are selected as seeds; it is never
expressed as a positive/negative sign on the seed's weight — see "seed
weighting" below for why that distinction is load-bearing, not cosmetic.

THE NETWORK AND ITS CONFIDENCE SCORE
-----------------------------------------
Same STRING interaction network as the Genomics module, filtered to
connections at or above a CONFIDENCE CUTOFF you choose (default 0.700, "high
confidence"). Only organisms STRING actually has a network for are offered;
picking one your uploaded genes don't otherwise support is refused up front
with an explanation, before any download.

STRING is identifier-agnostic (UniProt in Genomics, gene symbol/Ensembl/
Entrez here) internally: it keys every protein by its OWN "STRING protein
ID". This module bridges your identifiers onto that space itself — primarily
by gene SYMBOL (matched against STRING's own preferred name for each
protein), falling back to Ensembl gene ID and then Entrez ID for genes whose
STRING nomenclature differs from your table's (common in less-studied
organisms, where STRING sometimes uses locus tags instead of a common
symbol). The "Gene mapping" pane on the results screen reports exactly how
many genes matched via each tier, and refuses the run outright below 10%
overall mapping — a 2%-mapping run would otherwise "succeed" and return
noise with no indication anything was wrong.

STRING'S DATA SOURCE AND LICENCE
-------------------------------------
STRING (v12.0, as used by this app) is freely available under a Creative
Commons Attribution 4.0 (CC BY 4.0) licence, which requires giving STRING
appropriate credit and disclosing any changes made to the data. If you use
results from this module, cite STRING itself alongside Omiscera — see
string-db.org for the current recommended citation.

SEED WEIGHTING
-------------------
Two choices, both producing STRICTLY POSITIVE weights only:
  - EQUAL (default) — every seed weighs 1.0 (or a small integer multiple,
    if more than one row of your table mapped to the same protein — this
    happens with duplicate probes/transcripts for one gene, and their
    contributions are summed rather than kept separate). The padj/log2FC
    cutoffs already made the evidence call; re-weighting by the same
    statistics again would double-count it.
  - WEIGHTED BY |RANKING METRIC| — seeds are weighted by the ABSOLUTE value
    of whichever ranking metric GSEA uses (|log2FC|, |π-value|, or
    |-log10(p)|/|-log10(padj)|), winsorized at the 99th percentile among
    your seeds and floored, bounding the strongest seed's weight at exactly
    20x the weakest.

WHY ABSOLUTE VALUE, NEVER A SIGNED WEIGHT. The shared propagation engine's
seed-mapping step will shift any non-positive weight set upward by
|minimum| + 1 so every seed has strictly positive mass to spread (random-walk
seeds cannot be negative). Because the walk then normalizes seed mass to sum
to 1, that shift changes RATIOS, not just absolute values — a signed weight
would silently turn "all genes, both directions" into a run biased toward
whichever direction had the larger raw magnitude, with nothing on screen to
suggest it. Direction is therefore expressed ONLY by which genes are selected
as seeds (the "all/up/down" control above); every weight this module actually
produces is asserted strictly positive before propagation runs, so that shift
is designed to never fire here.

THE PERMUTATION NULL
-------------------------
Identical method to Genomics: a degree-preserving permutation null (each
seed is repeatedly swapped for a random protein of similar network degree,
many times), converted to a Z-score and a normal-tail p-value per protein,
then Benjamini-Hochberg FDR-corrected across the whole network. Permutation
count is adjustable (200/500/1000/2000; default 500) — deliberately capped
lower than Genomics' range, since every permutation here is a full walk over
the whole organism's network and this server runs several jobs concurrently.

A network below ~20 proteins at your chosen confidence is refused outright,
with an explanation: below that size, every "degree bin" the permutation
draws replacement seeds from holds only one candidate, so a permuted seed set
is IDENTICAL to the real one every time — the null could never move, and
every result would silently and uninformatively read as "not significant."

SEED CIRCULARITY — READ THIS BEFORE TRUSTING A "SIGNIFICANT SEED"
------------------------------------------------------------------------
A seed receives propagation signal DIRECTLY, every step, by construction. Its
null, meanwhile, is built from permutations in which it is usually NOT a
seed (only reselected occasionally, by chance, from its degree bin). That
comparison is structurally tilted in the seed's favour — a seed clearing the
significance threshold is therefore mostly EXPECTED, not a discovery.
Measured directly on this app's bundled human example (default settings,
500 permutations): of all seeds, roughly 97% clear FDR < 0.05 purely from
receiving their own signal, while under 0.1% of NON-seed proteins ever do —
and in one full run through the actual upload-to-results pipeline, 534 of the
571 significant proteins (93.5%) were the user's own seeds, and only 37 were
proteins never uploaded at all. Seed Z-scores can run into the hundreds.

The genuinely new information this module adds is the RECRUITED set —
significant proteins that were never seeds. The Hit table defaults to
showing RECRUITED genes only whenever a run produced any, specifically so
this small, actually-informative set is not buried under hundreds of
structurally-expected seed rows; a toggle switches to the full table (every
protein, seeds included), and the CSV download always contains everything
regardless of that toggle's state. This changes what the interface shows
first; it changes NOTHING about the underlying ranking, Z-scores or FDR.

MODULE DETECTION
---------------------
Same Markov Clustering (MCL) method as Genomics, run on the subnetwork
induced by the significant hits. Because a DEG list here is typically
hundreds to thousands of genes — not the tens Genomics seeds with — the hit
set is capped by Z-score (strongest signal first, default 600 genes) before
clustering, and the cap used is reported. Without this cap, an unconverged
MCL run on a large, dense subgraph tends to merge everything into one giant
"module" with a single generic label: entirely plausible-looking on screen
and entirely wrong. A "largest module contains most of the clustered genes"
warning fires for the same reason if it happens anyway.

Every module reports how many of its members were seeds versus recruited,
and — for members that were in your DE table — how many moved up versus
down. A module mixing both directions is not a single coherently-regulated
unit; read its label with that in mind.

MODULE LABELLING
---------------------
Each module is labelled by testing its own gene symbols against GO or
MSigDB (your choice; MSigDB covers human and mouse only). THE UNIVERSE FOR
THIS TEST IS THE NETWORK'S OWN ANNOTATED GENE SYMBOLS, not your DE table's
tested genes — a deliberate departure from the obvious default. Module
members are drawn from the network, so the network is the correct sampling
frame; using the DE table instead would silently delete every RECRUITED gene
from its own module's enrichment result — the exact failure this module
exists to avoid, and one that would look correct on screen while quietly
being wrong. Reactome is not offered here: Transcriptomics_module's Reactome
path builds its gene bridge from the DE table itself, which would reintroduce
exactly that same problem.

MULTIPLE-TESTING CORRECTION
-------------------------------
Propagation p-values are Benjamini-Hochberg FDR-corrected across every
protein in the network (default significance FDR < 0.05, adjustable without
re-running). Module labels are POST-SELECTION: modules are chosen for being
densely connected and carrying propagated signal, and only THEN tested for
enrichment — their FDR describes what a module contains, it is not a
hypothesis test, and it is not corrected for having compared multiple
annotation sources (GO's three namespaces, or several MSigDB collections)
per module.

RESULTS AND VISUALIZATION
-----------------------------
  - Overview — seeds, significant proteins, recruited proteins, modules, and
    network size at a glance, plus the measured null skewness (see below).
  - Gene mapping — coverage by identifier tier (symbol/Ensembl/Entrez),
    matching zero, one, or more than one protein (collisions, resolved by
    keeping the best-connected candidate and reporting the collision).
    Downloadable as CSV.
  - Hit table — every protein in the network, ranked by FDR: ID, gene
    symbol, log2FC/padj (blank for a recruited gene — it was never in your
    table), propagated score, Z-score, p-value, FDR, and whether it was a
    seed. Defaults to RECRUITED ONLY when the run produced any (see "seed
    circularity" above); the full table and the CSV download always contain
    every protein regardless. p-values/FDRs that genuinely underflow float64
    precision display as "< 1e-300" rather than a misleading "0.000e+0" —
    the CSV keeps full numeric precision either way.
  - Figures — the same interactive network and null-distribution figures as
    Genomics: node colour carries propagated score or module membership;
    node size only marks seeds apart from recruited proteins, it is not a
    measurement. Both exportable as SVG/PNG/TIFF.
  - Modules (when module detection ran) — one row per module: size, seed/
    recruited/up/down composition, mean log2FC among members with data, and
    top enriched term. Downloadable as CSV, alongside a module-coloured
    network figure.
  - Module enrichment — pick any module and see its full enrichment result
    (dot or bar plot, downloadable as its own CSV), not just the single top
    term used as the module's label.
  - Parameters used — every setting for this run, the single reference to
    check or screenshot for a methods section.

THINGS TO KEEP IN MIND
---------------------------
  - Seed circularity (above) is the most important caveat in this module: a
    "significant" seed is mostly expected, not discovered. Read the
    RECRUITED set as the finding.
  - p-values assume a normal null; the actual permutation null is measurably
    right-skewed (skewness around 2.0-2.5, measured on the bundled human
    example) — a caveat banner states the EXACT value for your run when it
    exceeds 1.0. Treat the RANKING as more trustworthy than the exact FDR.
  - Module labels are post-selection (see above) — a descriptive summary of
    what a module contains, not a hypothesis test of whether it should exist.
  - STRING's combined confidence score folds in text-mining and cross-
    organism evidence transfer alongside direct experimental evidence, so
    the network — and therefore this module's results — leans toward well-
    studied genes.
  - Less-studied, non-model organisms have both a sparser interaction network
    and sparser GO/MSigDB annotation, weakening propagation signal and module
    labels alike; watch the Gene mapping pane's coverage percentage.
  - A network below ~20 proteins at your chosen confidence is refused rather
    than silently returning an uninformative "nothing significant" result —
    try a lower confidence cutoff.
