GENOMICS — NETWORK PROPAGATION / MODULE DETECTION
=====================================================

WHAT THIS MODULE IS FOR
--------------------------
Network Propagation / Module Detection maps your genes onto the known
protein-protein interaction network (STRING) for your organism and spreads
("propagates") their signal across it, answering:

    "Starting from my genes as 'signal sources,' which other proteins does
    that signal reach and concentrate on as it spreads across the network
    — and do those proteins organize into densely-interconnected
    functional modules?"

This surfaces two kinds of results:
  - RECRUITED proteins — proteins you never uploaded, but that the network
    consistently pulls signal onto because they sit at a hub connecting
    several of your genes. Candidate "guilt by association" genes.
  - MODULES — dense neighborhoods of interacting proteins among your hits,
    which may or may not correspond to a single named pathway.

Two checkboxes control how far this goes:
  - "Network propagation" alone — run the signal-spreading step and report
    every protein whose accumulated signal is statistically significant.
  - "Module detection" — automatically includes propagation, and
    additionally clusters the hits into modules with a functional label
    for each. Results are additive: selecting Module Detection gives you
    everything Propagation gives you, plus the module table and figure.

USER INPUT
-------------
A table with:
  - One identifier column (required): Gene Symbol OR Ensembl Gene ID.
  - UniProt ID (optional) — skips resolution for that row if supplied.
  - A Gene Score (optional, numeric) — sets each seed's PROPAGATION
    STRENGTH (how much signal it injects), via a "Higher = more signal" /
    "Lower = more signal" toggle you choose alongside the column. If
    omitted, every uploaded gene is an equal-strength seed.

Genes are matched automatically once you select an organism and upload
your table; ambiguous matches are flagged for you to confirm. There's no
Homology (BLAST) option, since a raw sequence isn't a stable identifier to
anchor to a network node.

ORGANISM SELECTION
-----------------------
The organism list here is narrower than in other Genomics modules — only
species with an available protein-interaction network are offered, so
every organism you can pick is guaranteed to work.

CONFIDENCE CUTOFF
----------------------
The interaction network is filtered to connections at or above a
CONFIDENCE CUTOFF you choose (0.150 "low confidence" through 0.900
"highest confidence"; default 0.700 "high confidence," the standard
recommended default). A higher cutoff gives a smaller, cleaner, more
reliable network; a lower cutoff gives a larger network with more (but
noisier) connections. Genes with no counterpart in the network for your
organism are reported separately ("seeds not found in network") — they
simply can't participate.

STRING'S DATA SOURCE AND LICENCE
-------------------------------------
STRING (v12.0, as used by this app) is freely available under a Creative
Commons Attribution 4.0 (CC BY 4.0) licence, which requires giving STRING
appropriate credit and disclosing any changes made to the data. If you use
results from this module, cite STRING itself alongside Omiscera — see
string-db.org for the current recommended citation.

HOW PROPAGATION WORKS
--------------------------
Each of your genes starts as a "seed" holding some signal (1.0 by default,
or your Gene Score). That signal spreads across the network in a standard
"random walk with restart" process: at every step, signal on a protein
either moves to its network neighbors (weighted by interaction confidence)
or "teleports" back to the original seeds. This repeats until the signal
distribution stabilizes, producing a "propagated score" for every protein
in the network — not just your seeds.

WHY A PERMUTATION TEST, NOT A RAW SCORE CUTOFF
----------------------------------------------------
A raw propagated score isn't meaningful on its own: highly-connected "hub"
proteins accumulate a lot of signal from almost ANY set of seeds, simply
because they sit near many things — not because they're biologically
related to your genes. To tell these apart, the tool builds a null
distribution by re-running propagation many times with random seed sets
matched in network connectivity, then converts each protein's real score
into a Z-score and, from that null, a p-value, and finally corrects for
multiple testing (FDR) across the whole network. The number of
permutations is adjustable (500 to 10,000; default 1000). More
permutations give a more precise, stable result; 500 is the measured floor
below which the null's variance itself becomes unreliable to estimate — on
a real STRING-scale network, too few permutations can occasionally miss a
rare, high-leverage permutation and understate a protein's null spread,
manufacturing a spurious extreme Z-score rather than merely a noisy one.
2000+ is recommended for a final run or if borderline results look
unstable.

THE P-VALUE IS NOT A PLAIN NORMAL-TAIL APPROXIMATION. The permutation null
is measurably right-skewed (real human runs: skewness ~2-2.5), and a plain
normal approximation is wrong there by orders of magnitude in exactly the
tail where every hit is called. Each protein's p-value instead comes from a
Pearson III (shifted-gamma) fit to that protein's own null, by method of
moments, falling back to the normal tail only where the fit isn't
meaningful (little or no null variability, or too little measured skew to
fit reliably). An EMPIRICAL p-value (how many permutations actually met or
beat the observed score) is also reported for transparency, though it is
never itself used for FDR — at realistic permutation counts its floor
cannot survive correction across a whole network.

Proteins passing your significance threshold (default FDR < 0.05,
adjustable) are the "hits" — almost always including your original seeds,
plus any recruited proteins. A SEED clearing that threshold is largely
CIRCULAR: it receives propagated signal directly, while its null is built
from permutations in which it is usually not a seed, so the comparison is
structurally tilted in its favour. The Hit table therefore also reports a
second, stricter FDR computed within a family containing ONLY non-seed
proteins — the honest number for "is this recruited protein significant on
its own terms," not "significant in a family dominated by its own seeds."

MODULE DETECTION
---------------------
When selected, the significant hits are clustered into densely-
interconnected MODULES using Markov Clustering (MCL) — the same method
STRING itself uses. Modules smaller than your chosen minimum size (default
3 proteins) are dropped as too small to interpret.

Each module is automatically labeled by testing its own gene set against
Gene Ontology, giving it a human-readable "top enriched term" (e.g. "DNA
damage checkpoint") even though it was discovered purely from network
structure, not GO membership. Overly generic terms (annotated to more than
2000 genes) are excluded from labeling so the label stays specific and
useful.

The reported "Top enriched FDR" is CORRECTED for comparing multiple GO
namespaces (Biological Process, Molecular Function, Cellular Component) per
module: picking the smallest FDR out of several independently-tested
namespaces is itself a selection, and an uncorrected "best of N" reads more
significant than it is. The correction is min(raw FDR x namespaces actually
tested for that module, 1.0); the uncorrected, within-namespace value is
also shown for reference. Neither number is a hypothesis test — see
"Things to keep in mind" below.

RESULTS AND VISUALIZATION
-----------------------------
PROPAGATION (always shown):
  - Summary cards — seeds mapped, seeds not found in the network,
    significant hits, and recruited proteins. Check "seeds not found"
    first: if many of your genes never made it into the network, every
    result below is based on less input than you think.
  - Hit table — every protein tested, with its ID, propagated score,
    Z-score, p-value, an empirical p-value, FDR, a second RECRUITED-FAMILY
    FDR (blank for seeds — see above), and whether it was one of your
    seeds. Downloadable as CSV, full precision. Sort by FDR and look for
    rows where "Seed" is false: those are the recruited proteins, the
    genuinely new candidates this analysis exists to find — filter on the
    recruited-family FDR specifically for the honest version of that list.
  - Propagation network figure — the significant hit subnetwork, coloured
    by propagated score, with seeds outlined more heavily than recruited
    proteins. A recruited protein sitting in the middle of a dense cluster
    of your seeds is a stronger candidate than one hanging off the edge on
    a single edge.
  - Null-distribution histogram — every protein's Z-score, with a dashed
    line at the actual significance boundary. This shows you directly how
    far your hits sit above what random seed sets of the same degree would
    produce, rather than asking you to trust an FDR column on its own. If
    the line sits inside the bulk of the distribution rather than out in
    the tail, the result is weak no matter what the table says.

MODULE DETECTION (shown only when selected):
  - Module table — one row per module: size, members, top enriched term
    and its FDR. Downloadable as CSV.
  - Module network figure — the same subnetwork, recoloured by module
    membership instead of score. Use it to see whether your hits form one
    connected story or several unrelated sub-stories.
  - Per-module enrichment plots — pick any module and see its FULL
    functional enrichment as a dot or bar plot, not just the single top
    term used as the module's label. Worth opening whenever a module's
    label seems too generic to be useful: the second and third terms often
    describe the biology better than the first.

Both network figures are interactive (draggable nodes, adjustable font size
and node count, minimum edge confidence) and exportable as SVG, PNG, or
300 dpi TIFF. Only the significant hit subnetwork is ever drawn, so the
figure stays readable.

On the network figures, node COLOUR carries the meaning — propagated score
in the propagation figure, module membership in the module figure. Node SIZE
is not a measurement: it only marks seeds (larger, heavily outlined) apart
from recruited proteins (smaller). Don't read a bigger node as a stronger
result. On the per-module enrichment plots, circle size does encode gene
count, scaled so the largest circle is a fixed, readable size regardless of
how many genes are involved. Hover any node or point for its exact values.

A "Parameters used" panel on every results screen lists every setting
(both the ones you chose and the ones fixed by the tool) — the single
reference to check or screenshot for a methods section.

THINGS TO KEEP IN MIND
---------------------------
  - Seed circularity: a seed clearing the (all-protein) significance
    threshold is largely expected, not discovered — see "Why a permutation
    test" above. Read the recruited-family FDR, and the recruited proteins
    it flags, as the actual finding.
  - The network reflects functional-association evidence, not causality or
    direction — an edge means two proteins are associated, not that one
    regulates the other.
  - Recruited proteins are hypotheses, not confirmed hits — treat them as
    a "worth investigating" suggestion, not proven biology.
  - Hub bias is corrected for, but not perfectly — if results near the
    significance boundary look surprising, try more permutations.
  - Module labels are POST-SELECTION: modules are chosen for connectivity
    and signal, and only then tested for enrichment, so "Top enriched FDR"
    is not a hypothesis test even after the multi-namespace correction
    above.
  - Module detection only works on propagation's hit set — a real module
    that propagation didn't flag as significant won't appear. Loosening
    the confidence cutoff or FDR threshold gives module detection more to
    work with.
  - Per-module enrichment compares against the whole genome, not against
    your original gene list.
  - Best suited to hit-lists of at least a handful of genes that plausibly
    interact — a single isolated gene has no network to propagate through.
  - Less-studied, non-model organisms will have both a sparser network and
    sparser GO annotation, weakening both the propagation signal and the
    module labels.
