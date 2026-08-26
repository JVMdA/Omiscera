PROTEOMICS — QUANTITATIVE ENRICHMENT
=======================================

WHAT THIS MODULE IS FOR
--------------------------
Over-Representation Analysis (the other Proteomics enrichment module) starts
by drawing a line: proteins that passed your padj and fold-change cutoffs are
"significant", everything else is not. It then counts how many significant
proteins land in each pathway.

That works well when your signal is concentrated in a handful of strongly
changed proteins. It is blind, by construction, to the opposite situation:

    A whole complex or pathway shifts together — every member moves in the
    same direction by a modest amount — but no single protein moves far
    enough to clear the cutoff on its own.

To ORA, that pathway contains zero significant proteins, so it is never even
tested. Quantitative Enrichment answers a different question:

    "Are the proteins in this pathway systematically shifted — up or down —
    compared with every OTHER protein I quantified?"

Nothing is thresholded. Every quantified protein contributes its actual
measured value, and a pathway can come out significant because all of its
members moved a little, not only because a few moved a lot.

This is the method Perseus — the reference software in the proteomics field —
calls "1D annotation enrichment". Statistically it is a competitive
Wilcoxon rank-sum test, the same family as CAMERA and MAGMA, and the same
engine this app already uses for the Genomics Gene-set Burden Test.

WHY NOT GSEA?
-----------------
GSEA is the usual answer to "ORA throws away my magnitudes". It is not
offered here, deliberately.

GSEA assumes a long, dense ranked list — for RNA-seq that means ~20,000
genes, essentially the whole transcriptome. Shotgun proteomics typically
quantifies 3,000–8,000 proteins, and, more importantly, the ones it misses
are not missing at random: low-abundance proteins are systematically
under-detected. That biases one end of any ranking in a way GSEA's
permutation model does not account for.

The rank-sum test used here makes no assumption about list length or about
the shape of your value distribution, which is why the proteomics field
settled on it rather than on GSEA.

USER INPUT
-------------
The same resolved protein table the ORA module uses — you do not upload
anything twice. See README-ORA-proteomics.txt for the upload format,
protein-group handling (P04637;P04637-2), isoform stripping, and
contaminant/decoy filtering, all of which apply identically here.

The one difference: there is NO padj or fold-change cutoff on this screen.
Every quantified protein with a usable value is used.

THE RANKING METRIC
---------------------
You choose which measured quantity the test compares. Four options, the
same four the Transcriptomics GSEA module offers, so the definitions stay
consistent across the app:

  - log2 Fold-Change (DEFAULT) — the raw effect size. This is what Perseus
    1D annotation enrichment uses, and it is the right default for most
    proteomics questions: it answers "did this pathway's abundance shift?"
  - Signed -log10(p-value) — ranks by confidence, keeping direction.
  - π-value (log2FC × -log10 p) — combines size and confidence.
  - Signed -log10(padj) — like the second option, but using the adjusted
    p-value. Proteins with no padj are dropped and counted.

Proteins whose value cannot be computed (blank/non-numeric) are dropped and
reported to you as "Dropped (no finite score)". That count is per PROTEIN,
not per row — if one accession appears on several rows and any of them is
usable, that protein is scored, not counted as dropped.

If the same accession appears on more than one uploaded row (which happens
legitimately with protein groups), its values are AVERAGED rather than
letting whichever row came last silently win.

THE UNIVERSE (WHAT EACH PATHWAY IS COMPARED AGAINST)
--------------------------------------------------------
The comparison group is every OTHER protein you quantified — then narrowed,
per database, to those the database actually annotates. A protein mass
spectrometry never detected is not in the background, because it could never
have shown a shift in your experiment either.

This matters more in proteomics than transcriptomics, for the same reason
described at length in README-ORA-proteomics.txt: using a whole-proteome
background would treat thousands of never-measured proteins as if they had
been measured and found unchanged.

An optional uploaded background list narrows the universe further (e.g. to a
targeted panel) and applies to every selected database.

DATABASES
------------
The same eight sources as Proteomics ORA, so the two modules differ only in
the statistic, never in what they know:

  GO (Biological Process / Molecular Function / Cellular Component),
  Reactome, MSigDB, UniProt (disease), UniProt (keywords), InterPro,
  Complex Portal, DISEASES, and your own Custom database.

DISEASES is human-only and MSigDB is human/mouse-only; both grey out
automatically for other organisms rather than returning a silent empty
result. PANTHER is not included — see README-ORA-proteomics.txt for the
reason (an access barrier, not a licensing one).

THE STATISTICAL TEST
-----------------------
For each term, a TWO-SIDED Mann-Whitney U / Wilcoxon rank-sum test compares
the values of proteins inside the term against every other quantified,
annotated protein.

Two-sided is the important word, and it is the main way this module differs
from the Genomics Gene-set Burden Test that shares its engine. Burden scores
are unsigned severity values where only "higher" is meaningful, so that test
is one-sided. Proteomics log2FC is SIGNED: a complex being shut down is
exactly as interesting as one being switched on. A one-sided test would
silently miss every down-regulated pathway.

Each term reports how many of your proteins belong to it, the mean value
inside vs. outside (their difference is the "Effect"), a p-value, an
FDR-corrected p-value and a Bonferroni value.

The test is deterministic — no permutations, no random seed. The same input
always gives exactly the same answer.

READING THE "EFFECT" COLUMN
-------------------------------
Effect (mean difference) = mean inside the term − mean outside it.

Its SIGN is the direction, and this is where the biology is:
  - Positive → the term's proteins are higher than the rest of the
    quantified proteome (by your chosen metric).
  - Negative → they are lower.

This replaces ORA's "Fold Enrichment". They are not comparable numbers: fold
enrichment is a ratio of COUNTS ("3× more significant proteins than
expected"), while Effect is a difference in MEASURED VALUES ("this pathway
sits 1.4 log2 units below everything else").

Because direction lives in the Effect column, this module does NOT split
results into up / down / all the way ORA does. Each database produces ONE
result table, not three — which also means a materially smaller
multiple-testing burden.

MULTIPLE-TESTING CORRECTION
-------------------------------
Benjamini-Hochberg FDR is applied separately within each database (and each
selected MSigDB collection). There is NO correction across those families:
scanning every table and reporting the single best hit inflates your true
error rate above whatever any one FDR column says. A stricter Bonferroni
column is provided alongside. Significance defaults to FDR < 0.05 and is
adjustable on the results screen.

TERM-SIZE FILTERING
-----------------------
A minimum term size (default 3) excludes tiny, unstable sets; an optional
maximum (default 500) excludes very large, generic ones. A term also needs
at least one protein outside it to compare against — which is why a custom
database containing only ONE gene set cannot be tested: that set would be
the entire universe, leaving no comparison group.

RESULTS AND VISUALIZATION
-----------------------------
Each database produces a downloadable results table plus five figures.

1. PROTEIN SET OVERVIEW — "which pathways shifted, how far, and which way?"
  - Dot plot — one row per term. Horizontal position is the Effect, so
    terms left of zero shifted DOWN and terms right of zero shifted UP.
    Colour is significance; circle size is how many of your proteins the
    term contains.
  - Bar plot — the same terms as bars, rankable by FDR or by Effect.
    Ranking by Effect asks "where is the shift largest?"; by FDR, "where am
    I most confident?" These are different questions and often give a
    different top of the list.
  - Volcano plot — Effect against significance for EVERY tested term at
    once. Unlike the other four figures, this one draws all terms
    regardless of significance, so it stays informative even when nothing
    passes the threshold (see below).

2. PROTEIN SET RELATIONSHIPS — "are these separate findings, or one story?"
  - Enrichment map — links two significant terms when they share proteins.
    A tight cluster is usually ONE underlying signal reported under several
    names, not several independent discoveries. Needs at least two
    significant terms to draw anything.
  - Protein-concept network — shows the individual proteins behind each
    term, each coloured by its own log2FC. This is where you see whether a
    term is significant because many proteins shifted moderately or because
    a couple of extremes are dragging the set — a distinction the table
    cannot show you, and one that should change how much you trust the
    result.

All charts are interactive (draggable nodes, adjustable font size, number of
terms shown, and significance threshold) and exportable as SVG, PNG or
300 dpi TIFF. The results table downloads as CSV with full numeric
precision.

WHEN NOTHING IS SIGNIFICANT
-------------------------------
This is a normal outcome, not an error, and it is treated as one. If no term
passes your FDR threshold you still get the volcano plot — showing every
tested term and how close it came — plus a note saying how many terms were
tested. The other four figures are withheld, because they only draw
significant terms and would be empty.

If that happens, it is worth checking whether your pathways might be shifting
in BOTH directions at once (some members up, some down). Such a pathway has a
mean shift near zero and is genuinely invisible to this test by design — that
is precisely the case ORA is better at, and the two modules are complementary
for exactly this reason.

HOW THIS DIFFERS FROM ORA — AT A GLANCE
-------------------------------------------
                        ORA                    Quantitative Enrichment
  Question        "are significant hits    "did this pathway shift?"
                   over-represented?"
  Test            Fisher's exact           Wilcoxon rank-sum (two-sided)
  Uses            counts of significant    every quantified protein's
                  proteins                 measured value
  Cutoff          padj + |log2FC|          none
  Effect size     Fold Enrichment          Effect (mean difference)
  Direction       separate up/down/all     the SIGN of Effect
                  result tables
  Blind spot      coordinated but          pathways shifting in BOTH
                  sub-threshold shifts     directions at once

Running both on the same data is normal and informative. They can legitimately
disagree, and the disagreement itself is a result: a pathway that ORA finds
but this module does not is usually one whose members moved in opposite
directions.

THINGS TO KEEP IN MIND
---------------------------
  - Effect (a plain average difference) and FDR (from a rank-based test)
    usually agree but are not the same thing — read Effect as "how big" and
    FDR as "how confident".
  - The test treats each protein as independent. Members of the same complex
    are strongly co-regulated in reality, which can make such terms look
    slightly more significant than they truly are — a normal, standard
    trade-off for this class of fast test, shared with CAMERA and MAGMA.
  - Works best with a broad, representative set of quantified proteins,
    ideally everything your run could plausibly detect, rather than a
    hand-picked shortlist.
  - Overlapping database terms are not independent of each other. Half a
    dozen significant, heavily-overlapping GO terms usually represent one
    finding, not six — the enrichment map is there to make that visible.
  - A significant result means association, not causation.
