PROTEOMICS — PTM / PHOSPHOSITE ENRICHMENT
=====================================================

WHAT THIS MODULE IS FOR
--------------------------
Over-Representation Analysis, Quantitative Enrichment and PPI Network /
Module Detection all work on PROTEIN-level data - one row per accession.
This module is different: it works on SITE-level data - one row per
modified residue (e.g. "TP53 S315") - and asks a question none of the other
three can:

    "Given my phosphorylation changes, which KINASES look more or less
    active, based on how THEIR known substrate sites moved?"

    "Are my significantly changed sites over-represented among any
    kinase's known substrates?"

USER INPUT — SITE-LEVEL, NOT PROTEIN-LEVEL
-----------------------------------------------
This tile has its OWN upload, entirely separate from the other three tiles'
protein-level one - a phospho-enrichment analysis cannot start from a table
that has already collapsed every site on a protein down to one row.

Required: a protein identifier (accession, or a protein group like
"P04637;P04637-2") + EITHER a single combined site column ("S473") OR
separate residue + position columns (MaxQuant's own "Amino acid" +
"Positions within proteins") + log2FC + p-value. Optional: padj/FDR
(computed via Benjamini-Hochberg if omitted), multiplicity, localization
probability. MaxQuant's own "Phospho (STY)Sites.txt" export is accepted
directly.

ISOFORMS ARE NEVER STRIPPED. Every other Proteomics tile collapses
"P04637-2" to "P04637" because GO/Reactome/InterPro are annotated on the
canonical entry. This tile does NOT: "P04637-2" S100 is not necessarily the
same residue as "P04637" S100, so collapsing them would silently relabel a
site onto the wrong sequence. A separate "canonical accession" is tracked
alongside the real one, used ONLY to look up gene symbols and to match
against iPTMnet's own (canonical-only) kinase-substrate data - never to
merge two different rows' identities.

Ambiguous protein groups (more than one possible protein) are kept, not
dropped, but flagged - a site's position is only guaranteed valid for the
leading member of the group. Multiply-phosphorylated peptides (MaxQuant's
___1/___2/___3 classes at the same site) are collapsed to the
most-significant row, NEVER averaged - they are different analytes, and
averaging a doubly-phosphorylated peptide's quantification with a singly-
phosphorylated one at the same site would mix two different biological
measurements into one number. A localization-probability filter (Class I,
>= 0.75 by default) drops low-confidence site assignments when that column
is present; its absence is flagged, not silently ignored.

DATA SOURCE AND LICENCE
----------------------------
iPTMnet (research.bioinformatics.udel.edu/iptmnet), release 6.2 as used by
this app - the only source in this tile of kinase-to-substrate-SITE
relationships (~23,000 enzyme-substrate-site triples across 11 PTM types;
this tile uses only the Phosphorylation subset). Covers human, mouse, rat,
and a limited set of other model organisms - an organism iPTMnet has no
data for is refused with an explanation, not a silent empty result.

LICENCE: iPTMnet is Creative Commons Attribution-NonCommercial-ShareAlike
4.0 (CC BY-NC-SA 4.0) - stricter than every other data source this app uses
(STRING/GO/Reactome/DISEASES/Complex Portal are all CC BY 4.0 or CC0):
NON-COMMERCIAL USE ONLY, and any redistribution of iPTMnet's own data must
carry the same licence and attribution. This app's own source code ships
freely on GitHub without bundling any downloaded database (cache/ is never
committed - every fetcher re-downloads at each user's own install), so this
tile's code distribution itself is unaffected; but if you use results from
this tile, or you run this app for commercial purposes, you are personally
responsible for confirming your use complies with iPTMnet's licence. Cite
iPTMnet alongside Omiscera if you publish results from this tile.

KINASE ACTIVITY INFERENCE (ULM / KSEA-STYLE)
-------------------------------------------------
For each kinase, its known substrate sites' ranking-metric values (the same
metric choices as GSEA/Quantitative Enrichment - signed -log10(p) by
default, or log2FC, pi-value, signed -log10(padj)) are regressed against
every OTHER site you measured. A POSITIVE Activity Score means that
kinase's substrates moved UP (more phosphorylated) relative to the rest of
the experiment - the kinase appears MORE active; negative means less
active.

This reuses the exact same numerical engine as Transcriptomics'
Transcription Factor / Regulon "activity" analysis (a closed-form ULM
regression, a gene/site-label permutation null, and empirical-null
calibration) - see that tile's own README for the full statistical
derivation of each of the THREE reported columns (parametric / calibrated /
permutation), which applies here unchanged. ONE difference: unlike a TF,
which can ACTIVATE some targets and REPRESS others, every kinase-substrate
edge here is a phosphorylation event - there is no "this kinase represses
this site" concept in iPTMnet (a phosphatase is a completely different
enzyme, not a sign on the same edge). So there is no Activated/Repressed
split here; "Substrate Sites Tested" is reported instead.

IMPORTANT CAVEAT: known substrate counts are wildly uneven across kinases.
A handful (MAPK1/3, PRKACA, CSNK2A1, CDK1, AKT1, GSK3B, and a few others)
carry hundreds of curated substrate sites in iPTMnet; most kinases carry far
fewer. Statistical power is therefore NOT uniform across the results table,
and a top hit can reflect how thoroughly a kinase has been studied, not how
biologically important it is in your experiment. Read the table alongside
"Substrate Sites Tested" for each kinase, not the rank alone.

KINASE SUBSTRATE-SET ENRICHMENT
------------------------------------
A conventional over-representation test (Fisher's exact + Benjamini-
Hochberg FDR, the same engine every other enrichment tile in this app uses):
are your significantly changed sites (padj/log2FC cutoffs, tested
separately for all/up/down) over-represented among any kinase's known
substrate sites? Runs together with kinase activity in the same job, since
both use the same iPTMnet fetch.

MULTIPLE-TESTING CORRECTION
--------------------------------
Kinase activity: THREE separately-reported p-values per kinase (parametric,
permutation, calibrated), each independently BH-corrected across every
kinase tested - the calibrated FDR is the primary, sort-by column. Substrate
enrichment: standard BH-FDR within each direction's family, same as every
other ORA tile.

RESULTS AND VISUALIZATION
-----------------------------
  - Overview - sites/kinases tested, licence note.
  - Kinase activity - full sortable table (Activity Score, Calibrated Score,
    all three p-values/FDRs), downloadable as CSV, plus a diverging bar plot
    of the top kinases by |calibrated score|.
  - Kinase substrate enrichment - one table per direction (all/up/down),
    downloadable as CSV, with dot and bar plots (the same builders every
    other enrichment tile in this app uses).
  - Parameters used - every setting for this run.

THINGS TO KEEP IN MIND
---------------------------
  - Substrate-annotation bias (above) is the most important caveat here - a
    top-ranked kinase may simply be the best-studied one, not the most
    biologically active one in your data.
  - Isoform preservation (above) is deliberate and load-bearing - if you see
    the same protein with "-2"/"-3" suffixes carrying DIFFERENT sites,
    that is working as intended, not a data-quality problem to "fix" by
    merging them yourself before uploading.
  - This tile is currently the ONLY source in this app of kinase-substrate
    relationships. A future addition (dbPTM, not yet built) would add much
    deeper site coverage and a "known/novel site" flag, but has no clean
    kinase-substrate table of its own - it cannot replace iPTMnet here, only
    supplement it.
  - iPTMnet's CC BY-NC-SA 4.0 licence is more restrictive than this app's
    other data sources - see above before using results commercially.
