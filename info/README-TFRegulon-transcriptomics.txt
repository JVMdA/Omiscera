TRANSCRIPTOMICS — TRANSCRIPTION FACTOR / REGULON
====================================================

WHAT THIS MODULE IS FOR
--------------------------
A transcription factor (TF) is a protein that switches genes on or off by
binding their DNA. A "regulon" is simply the set of genes one TF controls.
This module takes your differential-expression results and asks two
different questions about the TFs behind them:

    "Is this TF's target set over-represented among my differentially
    expressed genes?" (TARGET ENRICHMENT)

    "Did this TF's targets move in the direction activating - or repressing
    - that TF would predict?" (ACTIVITY INFERENCE)

These are genuinely different questions. A TF can pass target enrichment
(its targets show up more than chance among your DE genes) without any clear
direction, and a TF can show strong activity (its targets moved in a
coherent, predicted direction) even if only a modest fraction of them were
individually "significant". Run either, or both.

USER INPUT
-------------
Shares the same differential-expression retrieval step as Over-Representation
Analysis and Gene Set Enrichment Analysis (see the Retriever README) - the
resolved table, with every gene's log2 fold-change, p-value and adjusted
p-value, is reused directly. No extra upload is needed for this module
itself; the regulon databases below ARE the gene-set source.

THE FOUR DATABASES
----------------------
  - TRRUST — curated TF-target interactions from the literature, each
    labeled Activation, Repression, or Unknown.
  - DoRothEA — TF-target interactions from OmniPath, each edge graded by
    confidence (A, strongest curated evidence, through E, weakest/purely
    computational). You choose how many grades to include via the "DoRothEA
    confidence levels" setting; A+B+C is the standard cut used by decoupleR
    itself, and loosening it gives broader coverage at the cost of noisier
    edges.
  - CollecTRI — DoRothEA's successor, from the same research group; broader
    coverage, also fetched via OmniPath and also SIGNED (each edge carries
    activation or repression). Note that the confidence-level control applies
    to DoRothEA only - CollecTRI is always used at its full coverage, so
    changing that setting will not alter your CollecTRI results.
  - ChEA3 — six ChIP-seq/co-expression-derived TF-target libraries. UNSIGNED:
    it records which genes a TF's binding is associated with, but not
    whether that binding activates or represses them.

TRRUST, DoRothEA and CollecTRI cover human and mouse. ChEA3's libraries use
human gene symbols and are human only; selecting a mouse organism disables
it, with an on-screen note explaining why - the same pattern used elsewhere
in Omiscera for a human-only database (e.g. MSigDB Hallmarks in ORA/GSEA).

Because ChEA3 carries no activation/repression information at all, it
participates in Target Enrichment only. Selecting it for Activity Inference
produces an explicit "no mode-of-regulation information" note instead of a
silently empty result.

TARGET ENRICHMENT
---------------------
Exactly the same statistical test as Over-Representation Analysis, pointed
at a TF's target set instead of a GO term or pathway: a one-sided Fisher's
exact test comparing how many of your differentially expressed genes belong
to the TF's targets against how many would be expected by chance, corrected
for multiple testing with Benjamini-Hochberg FDR (applied separately per
database). The padj/log2FC thresholds that define "differentially expressed"
here are the same two settings you would use in ORA. Direction of regulation
is ignored in this analysis - a TF's activated and repressed targets are
tested together as one set.

ACTIVITY INFERENCE
-----------------------
The one calculation specific to this module. Every tested gene gets a
continuous ranking value from the same metric choice GSEA uses (signed
-log10(p-value) by default). For each TF, that value is regressed against a
+1/-1/0 vector built from the TF's OWN signed targets: +1 for a gene the TF
is recorded as activating, -1 for one it represses, 0 for every gene that
isn't one of its targets. The resulting slope, standardized into a
t-statistic, is the TF's ACTIVITY SCORE (this is decoupleR's own "ulm"
method):

    - A clearly POSITIVE score means the TF's activated targets moved up and
      its repressed targets moved down relative to the rest of the
      experiment - i.e. the TF appears ACTIVE.
    - A clearly NEGATIVE score means the opposite pattern - the TF appears
      REPRESSED.

This is inferred from downstream mRNA movement, not measured directly - a
TF's own protein activity is often controlled after translation (by
phosphorylation, localization, degradation, etc.), so its mRNA level and its
real activity can disagree. The "TF activity vs. own expression" figure
(below) shows this directly: a TF read as active whose OWN mRNA barely moved
is a candidate for exactly this kind of post-translational activation.

A TF is only tested if enough of its targets - by count, not fraction - were
present in your experiment (the "Minimum regulon size" setting); below that,
a slope has too little support to trust and the TF is excluded rather than
reported.

THREE SCORES PER TF, AND WHY THERE ARE THREE
-------------------------------------------------
Every TF gets three views of its significance, computed three different ways:

  - PARAMETRIC — the ULM t-test itself (decoupleR's own method).
  - CALIBRATED — the parametric score, rescaled by a factor (λ, "lambda")
    measured from this run's own data (see below). THIS is what sorts and
    highlights the results table.
  - PERMUTATION — built by actually re-shuffling, thousands of times, which
    genes count as this TF's targets, and seeing how extreme a same-sized
    random set's slope typically is. Makes no assumption that the ranking
    values are normally distributed. Kept as a separate, deliberately
    UNcalibrated diagnostic.

WHY CALIBRATION EXISTS: a TF whose real targets are genuinely co-regulated
(move together for reasons beyond this one TF) can look more extreme than
the parametric t-test's null expects, because that test assumes every
gene's ranking value is an independent draw - co-regulated targets violate
that. Rather than leave this as an untested caveat, Omiscera measures it
directly, every run: since most TFs in a real dataset are not truly active,
the BULK of the observed activity scores approximates the true null - and
comparing that bulk's actual spread against plain statistical theory
reveals exactly how inflated the raw p-value is. That comparison IS λ:

    λ ≈ 1.0   the observed spread matches plain theory
    λ > 1.0   the observed spread is wider than theory by roughly that factor

λ is specific to YOUR data and to each database, and will differ from run to
run - there is no "normal" value to compare against. For a sense of scale
only: on this app's own bundled example dataset (a single simulated human
transcriptome), λ came out between roughly 1.4 and 1.7 across TRRUST,
DoRothEA and CollecTRI. That is one dataset, not a published expectation -
read your own λ, don't calibrate your intuition to that number.

READ THIS BEFORE TRUSTING λ - IT IS STILL AMBIGUOUS, LESS THAN IT WAS
--------------------------------------------------------------------------
A large λ has TWO possible causes, and this app cannot tell them apart with
certainty:

  (a) The null really is inflated - your TFs' targets are co-regulated, so
      the raw p-values were genuinely too optimistic. Calibration is right,
      and correcting is what you want.

  (b) A large fraction of your TFs are genuinely ACTIVE. Real signal widens
      the observed spread of scores too. Here there is nothing to correct,
      and calibration is more conservative than it needs to be.

This is measured, not hypothetical - and it drove a real change to how λ is
computed. A plain median/absolute-deviation estimate (the standard robust
choice) turned out to be fooled by (b) badly: a simulation with NO inflation
at all but 25% of TFs genuinely active produced λ ≈ 1.4-1.5 from that
estimator alone - indistinguishable from genuine inflation of the same size.
Omiscera now uses a TRUNCATED-NORMAL estimate instead: it looks only at the
central bulk of the score distribution, where real signal is least likely to
intrude, and corrects analytically for having discarded the tails. Measured
improvement, at the TF counts this module typically produces:

                              plain estimate   truncated estimate   (truth)
    0% truly-active TFs          ~1.0-1.1          ~1.0-1.1           1.0
    10%                          ~1.1-1.15         ~1.03-1.06         1.0
    25%                          ~1.4-1.5          ~1.04-1.1          1.0
    40%                          ~2.0              ~1.05-1.1          1.0

Two things had to hold for this to be a real fix rather than trading one
problem for another - both checked directly, not assumed:
  - it must still catch GENUINE inflation, or "fixing" the false alarm would
    destroy the feature: a true 1.6x inflation is still recovered at ~1.56,
    a true 3.0x at ~2.92.
  - under GRADED (continuous) signal, not neatly separated high/low TFs -
    the realistic shape - it is meaningfully better but not perfect: 40%
    graded activity gives λ ≈ 1.16 rather than 1.0, versus ≈ 1.43 with the
    plain estimator. The ambiguity is narrowed, not eliminated.

The calibration QQ plot does NOT resolve the remaining ambiguity - that was
tested directly. The idea was that real signal would leave the middle of the
plot on the y = x line while only the tails departed. It does not:
contamination in the tails compresses the rank axis and steepens the middle
too, so genuine inflation and real signal still look similar there. Telling
them apart with certainty would need the sample-level expression matrix this
app does not have - it works from a differential-expression summary table.

WHAT THIS MEANS FOR YOU. The calibrated column is the CONSERVATIVE default,
not the provably correct one. It is right under (a) and can still be a
little too strict under (b), though far less so than before this fix. It is
the default because a false positive in a published figure costs more than
a missed lead. But if you have independent reason to expect widespread
real TF activity - a strong perturbation, a well-powered experiment, a
pathway you already know is switched on - then read the PARAMETRIC column
instead and treat the calibrated one as a conservative lower bound. Both are
in the table, and in the downloadable CSV, for exactly this reason.

λ can only ever push results to be MORE conservative, never less - if the
observed spread came out narrower than theory, Omiscera does not "reward"
that by making results more significant.

BE PRECISE ABOUT WHAT THIS FIXES: calibration measurably REDUCES the excess
false-positive rate this problem causes, but does not eliminate it (in the
simulation used to validate this feature, the false-positive rate at a 5%
target dropped from about 19% to about 11% - closer to nominal, not exactly
nominal). It corrects for a run-wide inflation factor; it cannot separate a
single TF's own genuine co-regulation from every other TF's. Read a
calibrated result as "the best correction Omiscera can make from a
differential-expression summary table," not as a guarantee.

The permutation column is left deliberately UNcalibrated and carries the
SAME kind of blind spot as the parametric column, for a related but
distinct reason: it builds its null from RANDOM target sets, which are
typically uncorrelated, so it can under-represent how extreme a genuinely
co-regulated regulon looks. This is the same limitation already documented
for Gene Set Enrichment Analysis's own permutation test in this app (see the
GSEA README) - it is not specific to this module, and calibrating it too
would not have been justified without separately measuring it.

Practically: when calibrated and permutation agree, that agreement is
reassuring. When they disagree sharply for a TF you care about, that is a
sign worth investigating in the Regulon Detail figures below, not a cue to
pick whichever number is smaller.

MULTIPLE-TESTING CORRECTION
-------------------------------
Both analyses apply Benjamini-Hochberg FDR correction across every TF tested
in a given database, separately per database (activity inference corrects
the parametric, calibrated and permutation p-values separately, giving three
FDR columns). Default significance is FDR < 0.05, adjustable on the results
screen without re-running the analysis.

EACH RESULT TABLE IS ITS OWN FAMILY. Target enrichment splits by direction,
so selecting two databases already produces SIX independently-corrected
tables: TRRUST (all/up/down) and DoRothEA (all/up/down). Activity adds one
more family per database. There is no correction ACROSS those tables. So if
you scan all of them and report the best hit from whichever table it appeared
in, your effective false-discovery rate is higher than the FDR column
suggests - the number in the column is only valid within its own table. This
is the same convention Over-Representation Analysis uses (FDR per term
family, not global), and it is a deliberate choice, but it is easy to
forget when clicking between tables. Decide which table answers your
question BEFORE looking, rather than picking the table that gives the
prettiest result.

RESULTS AND VISUALIZATION
-----------------------------
TF ACTIVITY (one pane per selected signed database):
  - Summary cards — TFs tested, the ranking metric used, the number of
    permutations, the random seed, and λ ("Null inflation"). The seed makes
    a re-run with identical settings reproduce identical permutation and
    calibrated p-values exactly; quote the metric, permutation count, seed
    AND λ together if you report these results - λ is specific to this run's
    data and this database, and can differ run to run. If fewer than 30 TFs
    were tested, calibration is skipped (λ left at 1.0) and a note explains
    why - there is too little data in that database's run to estimate an
    empirical null reliably.
  - Calibration QQ plot — the figure that shows λ rather than asking you to
    trust the number. Every tested TF's activity score, sorted, plotted
    against the value a perfectly uninflated run would predict at that rank.
    A dotted "y = x" line marks "no inflation"; a solid line shows this run's
    ACTUAL fit. When the solid line visibly diverges from the dotted one,
    that divergence IS the inflation λ corrects for - worth opening before
    trusting any FDR column blindly, and useful even when nothing in this
    run reached significance, since it explains why.
  - Result table — TF, signed targets tested, how many of its targets are
    activated vs repressed, activity score, and three p-value/FDR sets
    (calibrated, permutation, parametric). Downloadable as CSV.
    "Signed targets tested" is the number of that TF's targets that are BOTH
    present in your data AND carry a direction - i.e. the genes the activity
    score was actually computed from. It always equals activated + repressed.
    Note it is NOT the same thing as the enrichment table's K, which counts
    the TF's whole target set in the universe; the two tables sit next to
    each other, so the columns are named differently on purpose.
  - TF activity bar plot — one horizontal bar per significant TF (by
    calibrated FDR), sorted by score. Red bars are TFs read as active, blue
    bars as repressed - the module's headline figure, meant to be read at a
    glance: which programs are switched on, which are switched off.
  - TF activity volcano — activity score (how strongly active/repressed)
    against -log10(calibrated FDR) (how confident, after correcting for this
    run's own measured inflation). Point size scales with the number of
    signed targets tested, so a
    lone extreme point that is also small is immediately visible as resting
    on very few targets - worth a second look before trusting it.
  - TF activity vs. own expression — activity score (x) against the TF's OWN
    log2 fold-change (y) in your data. Three things are encoded, and the
    figure carries this same key on screen:
      - FILLED marker = the TF's ACTIVITY is significant (calibrated FDR
        below your threshold); red = inferred active, blue = inferred
        repressed. Hollow grey = not significant.
      - DARK RING = the TF's OWN mRNA change is significant (its own padj
        passes the threshold). No ring = its own expression did not
        significantly change. A TF whose own padj is missing gets no ring -
        absence of evidence is not drawn as confidence.
      - SIZE = how many of that TF's targets are present in your data.
    The case worth hunting for is a FILLED RED circle with NO RING sitting
    near the horizontal zero line: the TF appears to be acting, but its own
    transcript never moved. That is the signature of activation happening
    after translation (phosphorylation, localisation, degradation) rather
    than by making more of the transcription factor - the caveat stated in
    words above, made visible.

REGULON DETAIL (one TF at a time, for the TFs with enough evidence to be
worth this level of detail - see "Minimum regulon size" above):
  - ULM regression plot — the actual regression behind that TF's activity
    score, drawn. Each target is plotted at its weight (+1 activated, -1
    repressed) against ITS OWN RANKING-METRIC VALUE - the exact number the
    fit consumed, which is why the y-axis is labelled with whichever ranking
    metric you chose rather than with log2 fold-change. The grey violin at
    weight 0 is every tested gene: it is most of the data, and it is what
    pins the line's height, so it is drawn rather than left implicit.
    Behind the fitted line, a fan of faint grey lines shows the SAME fit
    recomputed on random relabelings of which genes count as this TF's
    targets - a visual companion to the permutation p-value above, not a
    stronger proof than it (see "two p-values" above: the fan shares that
    method's own blind spot for genuinely co-regulated targets).

    One thing to expect: the fitted line does NOT have to pass through the
    middle of the activated or repressed points, and for a lopsided regulon
    (say 1 activated target and 6 repressed) it can pass nowhere near one of
    them. That is correct, not a drawing error - with ~99% of genes sitting
    at weight 0, that central mass pins the line's height and the few
    targets only tilt it. Reading the line together with the violin is the
    point of the figure.
  - Coherence bars — of this TF's activated targets, what fraction actually
    moved up; of its repressed targets, what fraction actually moved down.
    Separates COHERENCE from MAGNITUDE: a TF can have a large activity score
    driven by a few big movers, or a smaller score with near-unanimous
    agreement across many targets - the result table alone cannot show you
    which.

TF TARGET ENRICHMENT (one pane per selected database x direction):
  - Result table — the same columns as ORA (term/gene counts, fold
    enrichment, p-value, FDR, overlapping genes). Downloadable as CSV.
  - Dot plot and bar plot — the identical builders ORA uses elsewhere in
    Omiscera, so they read the same way: circle size (dot plot) is scaled so
    the largest circle in a chart is a fixed, readable size regardless of
    how many genes you uploaded, and AREA - not diameter - tracks the
    overlap count.
  - TF-TF overlap map — draws a line between two TFs when they share
    DIFFERENTIALLY EXPRESSED targets (the genes in the "Overlapping Genes"
    column, not their full recorded regulons). TF regulons overlap heavily
    in practice (two well-studied TFs can easily share a hundred targets) -
    a tight cluster of TFs in this map is usually ONE regulatory programme
    reported under several names, not several independent findings. The
    "Edge threshold" control sets how many shared DE targets are needed
    before a line is drawn.
  - TF-target network — shows the individual target genes behind each
    enriched TF, coloured by their own log2 fold-change, and which targets
    are shared between TFs.

A "Rank & show by" control on the enrichment charts lets you sort/display by
FDR, fold enrichment, gene count, or gene ratio - the same convention as
every other Omiscera dot/bar plot.

THINGS TO KEEP IN MIND
---------------------------
  - This module is effectively human-and-mouse only (ChEA3: human only) -
    unlike most of Omiscera, which works for any organism with adequate
    database coverage.
  - A regulon reflects literature and ChIP-seq evidence pooled across many
    studies and cell types, not your specific tissue or condition - a TF's
    real targets in your experiment may be a subset of, or differ from,
    its recorded regulon.
  - Calibration reduces, but does not eliminate, the co-regulation-driven
    optimism in the activity p-value - see "Three scores per TF" above. The
    permutation column is a separate, deliberately uncalibrated check, not a
    fallback "more correct" answer. All three are reported so you can see
    where they agree and where they don't, not so you can pick whichever
    number is smallest.
  - Target enrichment and activity inference can legitimately disagree: a TF
    can be enriched without a clear direction, or show a directional signal
    without every individual target reaching significance on its own.
  - A significant activity score is a hypothesis about which regulatory
    program is active, not proof that TF itself is the cause - as with any
    association-based result, treat it as a lead for follow-up, not a
    conclusion.
  - Double-check "Minimum regulon size" if a TF you expected to see is
    missing - a real, biologically important TF with few targets detectable
    in your experiment is excluded by design, not by mistake.
