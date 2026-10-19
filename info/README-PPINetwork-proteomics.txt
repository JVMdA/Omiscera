PROTEOMICS — PPI NETWORK / MODULE DETECTION
=====================================================

WHAT THIS MODULE IS FOR
--------------------------
Over-Representation Analysis and Quantitative Enrichment both ask "which
pre-defined gene set is over-represented / concentrated at the top of my
ranking?" This module asks two different questions, using the known
protein-protein interaction network (STRING) for your organism instead of a
gene set:

    "Starting from my significantly changed proteins as 'seeds,' which
    proteins does their signal concentrate on as it spreads across the
    network — INCLUDING proteins I never uploaded?"

    "Do the significant proteins organize into densely-interconnected
    functional MODULES?"

The numerics behind this — the network parser, the signal-spreading step,
the statistical test, and the clustering — are the SAME engine the Genomics
"Network Propagation / Module Detection" and Transcriptomics "Network /
Module Detection" modules use (see either module's own README for a
from-scratch walkthrough of the maths). Only one copy of that engine exists
(Genomics_network.py), and the orchestration around it
(Transcriptomics_network.run_network_analysis) is shared too — this module
supplies only the UniProt-accession-space identifier bridge and module
labelling, as optional hooks into that same orchestrator, so this tile and
Transcriptomics' Network tile can never quietly drift into disagreeing about
how propagation, the permutation null, or clustering itself work.
Everything documented below is what is DIFFERENT here:

    Transcriptomics                  Proteomics
    -------------------------------  ------------------------------------
    gene symbol / Ensembl / Entrez   UniProt accession (+ gene symbol fallback)
    tx.select_degs                   Proteomics_module.select_significant_proteins
    module labels: GO / MSigDB       module labels: GO / Reactome / Complex Portal

USER INPUT
-------------
Shares the same protein-level retrieval step as Over-Representation Analysis
and Quantitative Enrichment (see the Retriever README) — the resolved table
(UniProt accession + log2FC + p/padj) is reused directly. No extra upload is
needed.

Seeds are chosen from that table using the SAME padj/log2FC cutoffs and
"all/up/down" direction control as ORA and Quantitative Enrichment
(Proteomics_module.select_significant_proteins — the exact same call, so
this module can never disagree with the other two tiles about what counts as
a significantly changed protein). A DIRECTION choice restricts WHICH proteins
are selected as seeds; it is never expressed as a positive/negative sign on
the seed's weight — see "seed weighting" below for why that distinction is
load-bearing, not cosmetic.

THE NETWORK AND ITS CONFIDENCE SCORE
-----------------------------------------
Same STRING interaction network as Genomics and Transcriptomics, filtered to
connections at or above a CONFIDENCE CUTOFF you choose (default 0.700, "high
confidence"). Only organisms STRING actually has a network for are offered;
picking one your uploaded proteins don't otherwise support is refused up
front with an explanation naming the organism, before any download.

THE IDENTIFIER BRIDGE — WHY IT IS ITS OWN CODE PATH
---------------------------------------------------------
STRING keys every protein by its OWN "STRING protein ID" internally. This
module bridges your UniProt accessions onto that space via STRING's own
UniProt_AC alias file — the SAME kind of alias STRING publishes for gene
symbol/Ensembl/Entrez, just a different column — falling back to gene symbol
for the (typically small) remainder. Measured on the bundled human example:
about 95% of uploaded accessions map via the accession tier directly.

Two collision hazards apply to UniProt aliases specifically and are resolved
explicitly, never by a silent "last one wins" overwrite (the pattern
Genomics_network's own convenience index uses internally, which would
silently send a seed — or a module's enrichment evidence — to an arbitrary
protein): measured on the cached human STRING v12 alias file, 230 accessions
map to more than one STRING protein, and 17,224 STRING proteins carry more
than one UniProt accession (up to 267 for one heavily cross-referenced
protein). Resolution order, identical to the rule Transcriptomics' own
identifier bridge already uses for gene symbol/Ensembl/Entrez: only proteins
actually present in the graph are candidates; the highest-weighted-degree
protein wins (the better-connected protein is the one the network can
actually propagate through); ties break on the lexicographically smallest
STRING id, so the result is reproducible across runs. Every collision is
reported, never silently dropped.

The "Protein mapping" pane on the results screen reports exactly how many
proteins matched via each tier (accession vs. gene symbol), and refuses the
run outright below 10% overall mapping — a 2%-mapping run would otherwise
"succeed" and return noise with no indication anything was wrong.

STRING'S DATA SOURCE AND LICENCE
-------------------------------------
STRING (v12.0, as used by this app) is freely available under a Creative
Commons Attribution 4.0 (CC BY 4.0) licence, which requires giving STRING
appropriate credit and disclosing any changes made to the data. If you use
results from this module, cite STRING itself alongside Omiscera — see
string-db.org for the current recommended citation.

SEED WEIGHTING
-------------------
Two choices, both producing STRICTLY POSITIVE weights only — identical
scheme to Transcriptomics' Network tile (see that module's README for the
full derivation of why a signed weight would silently bias the whole run):
  - EQUAL (default) — every seed weighs EXACTLY 1.0, regardless of how many
    rows of your table happened to map onto the same STRING protein.
  - WEIGHTED BY |RANKING METRIC| — seeds are weighted by the ABSOLUTE value
    of whichever ranking metric Quantitative Enrichment uses (|log2FC|,
    |π-value|, or |-log10(p)|/|-log10(padj)|), winsorized and floored,
    bounding the strongest seed's weight at exactly 20x the weakest. A
    protein backed by more than one row takes the MAX magnitude among them,
    not the sum.

Direction is expressed ONLY by which proteins are selected as seeds (the
"all/up/down" control above); every weight this module actually produces is
asserted strictly positive before propagation runs.

THE PERMUTATION NULL
-------------------------
Identical method to Genomics and Transcriptomics: a degree-preserving
permutation null, with each protein's p-value from a Pearson III (shifted-
gamma) fit to its own null where the fit is meaningful, falling back to the
normal tail otherwise — see the Genomics module's README for the full
derivation; all three modules share one engine.

Permutation count is adjustable (500/1000/2000; default 1000) — see the
Transcriptomics Network README for the calibration evidence behind the 500
floor (a lower count was measured to manufacture spurious extreme Z-scores
from a single missed high-leverage permutation, not merely add noise).

A network below ~20 proteins at your chosen confidence is refused outright,
for the same reason: below that size, every "degree bin" the permutation
draws replacement seeds from holds only one candidate, so the null could
never move.

Propagation p-values are Benjamini-Hochberg FDR-corrected TWICE: once across
every protein in the network (the "FDR" column), and independently a SECOND
time across only the non-seed proteins (the "FDR (recruited family)"
column) — see "seed circularity" below.

SEED CIRCULARITY — READ THIS BEFORE TRUSTING A "SIGNIFICANT SEED"
------------------------------------------------------------------------
A seed receives propagation signal DIRECTLY, every step, by construction, so
a seed clearing the significance threshold is mostly EXPECTED, not a
discovery — see the Transcriptomics Network README's "seed circularity"
section for the full explanation (identical mechanism here; only the
identifier space differs). The genuinely new information this module adds
is the RECRUITED set — significant proteins that were never seeds, judged by
their OWN, seed-free FDR family. The Hit table defaults to showing RECRUITED
proteins only whenever a run produced any; a toggle switches to the full
table, and the CSV download always contains everything regardless.

A properly-calibrated result CAN legitimately be zero hits — that is not a
bug, it means the differential-abundance signal is not concentrated enough
on the STRING network to survive honest correction.

MODULE DETECTION
---------------------
Same Markov Clustering (MCL) method as Genomics and Transcriptomics, run on
the subnetwork induced by the significant hits, capped by Z-score (strongest
signal first, default 600 proteins) before clustering — see the
Transcriptomics Network README for why the cap exists and what the
"largest module contains most of the clustered proteins" warning means.

Every module reports how many of its members were seeds versus recruited,
and — for members that were in your table — how many moved up versus down.

MODULE LABELLING — AND WHY REACTOME WORKS HERE
---------------------------------------------------
Each module is labelled by testing its own STRING proteins against GO,
Reactome, or Complex Portal (your choice; default GO). THE UNIVERSE FOR THIS
TEST IS THE NETWORK'S OWN ANNOTATED PROTEINS, not your uploaded table's
tested proteins — the same rule Transcriptomics' Network tile uses, for the
same reason: module members are drawn from the network, so the network is
the correct sampling frame, and using the uploaded table instead would
silently delete every RECRUITED protein from its own module's enrichment
result.

Each database's own gene set (keyed by UniProt accession, since that is how
GO/Reactome/Complex Portal publish their data) is translated into
STRING-protein-ID space through the SAME accession bridge described above,
BEFORE testing — so a protein carrying several accessions (see the collision
numbers above) is counted exactly once per term, never once per accession.
Testing directly in accession space, the way the Genomics module's own ORA
does, would double-count exactly the proteins the collision numbers above
show are common.

REACTOME AND COMPLEX PORTAL ARE OFFERED HERE, UNLIKE IN TRANSCRIPTOMICS'
NETWORK TILE: Transcriptomics' identifier bridge is built FROM the uploaded
differential-expression table itself, so a recruited gene — one the bridge
never saw, by definition — cannot be found in Reactome's own gene sets
either, and Reactome is withheld there rather than silently under-counting.
This module's bridge is anchored on STRING's own UniProt accessions,
independent of what you uploaded, so a recruited protein's accession is
still resolvable and Reactome (and Complex Portal, which has no
Transcriptomics-side equivalent at all) work correctly. Reactome additionally
requires a resolved organism NAME (not just a taxon id) to filter pathways
by species; retrieval guarantees this is available whenever an organism was
selected.

MULTIPLE-TESTING CORRECTION
-------------------------------
Identical scheme to Transcriptomics' Network tile: propagation p-values are
Benjamini-Hochberg FDR-corrected across every protein in the network, and
independently again across the non-seed proteins only (default significance
FDR < 0.05, adjustable without re-running). Module labels are
POST-SELECTION — modules are chosen for being densely connected and carrying
propagated signal, and only THEN tested — so the "Top enriched FDR" shown is
the best within-namespace FDR multiplied by how many namespaces were tested
for that module (capped at 1.0); the raw, single-namespace value is also
shown as "Top enriched FDR (raw, uncorrected)" for comparison.

RESULTS AND VISUALIZATION
-----------------------------
  - Overview — seeds, significant proteins, recruited proteins, modules, and
    network size at a glance, plus the measured null skewness.
  - Protein mapping — coverage by identifier tier (accession, then gene
    symbol), matching zero, one, or more than one protein (collisions,
    resolved by keeping the best-connected candidate and reporting the
    collision). Downloadable as CSV.
  - Hit table — every protein in the network, ranked by FDR: STRING ID, gene
    symbol, log2FC/padj (blank for a recruited protein — it was never in
    your table), propagated score, Z-score, p-value, empirical p-value, FDR,
    FDR (recruited family), whether it was a seed, and seed score. Defaults
    to RECRUITED ONLY when the run produced any; the full table and the CSV
    download always contain every protein regardless.
  - Figures — the same interactive network and null-distribution figures as
    Genomics/Transcriptomics: node colour carries propagated score or module
    membership; node size only marks seeds apart from recruited proteins.
    Both exportable as SVG/PNG/TIFF.
  - Modules (when module detection ran) — one row per module: size, seed/
    recruited/up/down composition, mean log2FC among members with data, and
    top enriched term with both its corrected and raw FDR. Downloadable as
    CSV, alongside a module-coloured network figure.
  - Module enrichment — pick any module and see its full enrichment result
    (dot or bar plot, downloadable as its own CSV), not just the single top
    term used as the module's label.
  - Parameters used — every setting for this run, the single reference to
    check or screenshot for a methods section.

THINGS TO KEEP IN MIND
---------------------------
  - Seed circularity (above) is the most important caveat in this module: a
    "significant" seed is mostly expected, not discovered. Read the
    RECRUITED-FAMILY set as the finding, not the all-protein FDR.
  - A properly-calibrated run CAN legitimately find zero significant
    proteins — that means the signal was not concentrated enough on the
    network to survive honest correction, not that something is broken.
  - p-values are NOT a plain normal-tail approximation — see "the
    permutation null" above; treat the RANKING as more trustworthy than any
    single exact FDR.
  - Module labels are post-selection — a descriptive summary of what a
    module contains, not a hypothesis test of whether it should exist — but
    the namespace-selection bias from comparing multiple annotation sources
    per module IS corrected in the headline FDR shown.
  - STRING's combined confidence score folds in text-mining and cross-
    organism evidence transfer alongside direct experimental evidence, so
    the network — and therefore this module's results — leans toward well-
    studied proteins.
  - Less-studied, non-model organisms have both a sparser interaction network
    and sparser GO/Reactome/Complex Portal annotation, weakening propagation
    signal and module labels alike; watch the Protein mapping pane's
    coverage percentage.
  - A network below ~20 proteins at your chosen confidence is refused rather
    than silently returning an uninformative "nothing significant" result —
    try a lower confidence cutoff.
