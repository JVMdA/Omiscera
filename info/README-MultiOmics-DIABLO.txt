Multi-omics — DIABLO / Multiblock Integration
==============================================

What this method is for
------------------------
The other Multi-omics methods are UNSUPERVISED or descriptive. DIABLO (Data Integration Analysis for Biomarker
discovery using Latent cOmponents; Singh et al. 2019, Bioinformatics 35:3055 — mixOmics' block.splsda) is
SUPERVISED: given a GROUP for every sample (control / treated, two or more), it finds a few COMPONENTS per omics
that are:
  - correlated with one another ACROSS the omics (the omics "agree" on these components), and
  - discriminant: together they separate the groups,
while keeping only a SMALL, chosen set of features per omics (sparse). It is the supervised counterpart of MOFA:
where MOFA asks "what varies together" with no groups involved, DIABLO asks "what, together, tells the groups
apart".

THE CENTRAL PITFALL, AND WHY THE RESULT IS A CROSS-VALIDATED NUMBER
----------------------------------------------------------------------
A sparse model fitted and judged on the SAME samples it was optimised on will always look like it separates the
groups — that is what fitting does, not evidence of a real signal. The only honest number is how well the model
classifies samples it did not see when its features were chosen: repeated, stratified cross-validation picks the
number of components and how many features to keep per omics, and reports the cross-validated BALANCED ERROR RATE
(BER — the mean, over the groups, of each group's own misclassification rate; robust to unequal group sizes,
unlike plain accuracy). On top of that, Omiscera refits and re-cross-validates several times with the group
labels SHUFFLED (same data, same settings) and reports the share of those shuffles that did at least as well as
the real labels — a permutation p-value for "could an error this low happen on this data by chance". Measured
directly while building this: averaged over many INDEPENDENT random datasets the cross-validated error is unbiased
at the 50% chance level (two balanced groups), but on a SINGLE FIXED dataset one label permutation can look
deceptively good (varied from about 30% to 70% error in testing) simply from that dataset's own feature
correlation structure — which is exactly why a permutation test, not a single cross-validated number, is the
honest check.

WHAT IS HERE AND WHAT IS NOT
------------------------------
DIABLO's only published implementation is the R package mixOmics (Rohart et al. 2017, PLoS Comput Biol 13:
e1005752). Its algorithm — sparse generalized canonical correlation analysis (sGCCA; Tenenhaus et al. 2014,
Biostatistics 15:569), with the group dummy as an extra "block" in regression mode — is implemented here in pure
Python, written from mixOmics' OWN SOURCE, not from the paper alone (the paper does not specify every numerical
choice: the SVD initialisation, which block's score initialises the outcome's weights, the exact soft-thresholding
rule, which score deflates which block, the renormalisation used at a missing cell, or the three classification
distances). Every one of those choices was checked against a real local mixOmics run: fitting, prediction, all
three classification distances and both vote rules agree with mixOmics to float64 precision (differences of 1e-14
to 1e-15) — see tests/test_multiomics_diablo.py and tests/r_reference/README.txt. No R is needed to RUN Omiscera;
only to regenerate that check.

How to use it
--------------
1. Open the Multi-omics box "DIABLO / Multiblock Integration". It opens the same upload screen as Correlation and
   MOFA: choose two or three of Transcriptomics, Proteomics, Metabolomics (Genomics cannot take part — a gene list
   has no values per sample), pick the organism, upload one table per omics, and "Resolve identifiers and match
   samples". Formats and sample matching are exactly those of Correlation-based Integration (see
   README-MultiOmics-Correlation.txt). The "Example file" links give the same synthetic 36-sample example MOFA
   uses (Control / Treated, a planted treatment signal shared by all three omics, plus a cholesterol and a lipid
   signal that are real but unrelated to the groups — DIABLO should select the treatment features, not those).
2. "Continue to DIABLO setup": choose the sample-group source (REQUIRED — a Group column found in a file, or a
   sample sheet you upload), the model settings, and the feature-count grid to search. "Run DIABLO".
3. Results panes: Overview, Components, Selected features, Performance, Pathway enrichment (run on request, per
   omics), Data used, Parameters used. Every figure has the zoom / pan / in / out / autoscale / reset buttons, a
   Size menu, and "Download figure" (SVG / PNG / TIFF); every table can be downloaded as CSV.

Sample size
------------
DIABLO needs at least two groups with 3 or more samples each. Cross-validation needs more: below about 20 samples
total the cross-validated error is noisy and the selected features unstable — the setup page and the results say
so. A cross-validation fold cannot be larger than the smallest group (asking for more folds than that is refused).

Data preparation
------------------
For every omics: features missing in more than half of the samples are dropped; constant features are dropped;
the N most variable features are kept (default 1000 — DIABLO's authors advise pre-filtering to the most variable
features before fitting, same as MOFA); optional log2(x + 1) for raw counts / intensities. Inside the model itself,
every feature is centred and scaled to unit variance (mixOmics' own scale = TRUE), refitted on each cross-
validation training fold — never on the whole dataset, which would leak information from the held-out samples.

The model
----------
For each component, DIABLO looks for a weighted sum of a few features per omics (the "loading") such that the
resulting per-sample score is, as the DESIGN says it should be, correlated with the OTHER omics' scores and with
the group — sparse, because only `keepX` features per omics and component are kept (the rest are set to exactly
zero by soft-thresholding; the model is genuinely using only those features, not merely ranking all of them).

DESIGN: how strongly every PAIR of omics is linked while fitting. The outcome (the groups) is always linked to
every omics at 1 (mixOmics' own convention, confirmed from its source); the omics-to-omics link is a setting
(default 0.1, the DIABLO paper's own case-study choice) — 0 drops omics agreement and tends to pick features
specific to a single omics; 1 fully couples the omics and can pick features correlated with each other rather than
with the group.

Choosing the number of components and keepX. Neither is typed in blind:
  - the NUMBER OF COMPONENTS is chosen by cross-validated BER on a NON-SPARSE fit (every feature kept — the
    question here is only how many components help, not which features), one-standard-error rule: the smallest
    number of components whose BER is not worse than the best one by more than its own standard deviation;
  - keepX is then tuned ONE COMPONENT AT A TIME (mixOmics' own convention, and far cheaper than searching every
    component's grid jointly): with the keepX of earlier components fixed at their chosen values, every
    combination of the offered feature counts across the omics is cross-validated for that one component, and the
    combination with the lowest BER is kept.

Classifying a sample. Each omics' component gives a new sample a score; the predicted class comes from a DISTANCE
in score space (mixOmics' own three rules, chosen in the setup page, all computed and the chosen one shown):
  max.dist          the highest value of a per-component linear regression of the group dummy.
  centroids.dist     nearest class-mean score, Euclidean, using the components cumulatively.
  mahalanobis.dist   the same, scaled by the within-class score covariance.
The omics then VOTE on the class: majority vote, or a vote weighted by how well each omics' components tracked the
outcome's own fitted component on the training data (mixOmics' get.weights) — the weighted vote is what the
cross-validated error and the permutation test use.

What the selected features tell you
--------------------------------------
A selected feature is one the sparse model found USEFUL for telling the groups apart, together with the features
of the other omics — not necessarily THE biologically important one. Correlated features compete for a slot: a
feature can be left out only because another, similarly informative one happened to be picked instead. The design
value, the number of components and the folds are all choices; different reasonable choices can select different
features. Read a single run's selection as a lead to test, not a final list — the permutation p-value answers
"could labels this good happen by chance on this data", not "is this biology correct".

Pathway enrichment of the selected features
----------------------------------------------
Per omics, on request: are the SELECTED features over-represented in a Reactome pathway, against every feature
that entered the model (after the most-variable filter) as the background? One-sided Fisher exact test — the same
convention the single-omics ORA tiles use (the app's own per-pathway function, not re-implemented), Benjamini-
Hochberg across the pathways tested. Features are matched to a pathway by gene (transcripts and proteins, via
Entrez Gene ID / UniProt accession) or by compound name; features that match nothing are left out of the test. Run
per omics, so the pathway sets (and, for metabolites, the compound-name matching) are built only when asked — the
first metabolite request can take a few minutes, later ones reuse it.

What the numbers do not tell you
-----------------------------------
  - A component is a statistical pattern fitted to discriminate your groups, not a pathway and not a cause.
  - The cross-validated BER has its own uncertainty (shown as +/- a standard deviation over the repeats) — with
    few samples, do not read it as a single clean number.
  - The permutation p-value is about THIS dataset's labels, not a claim about the underlying biology.
  - Variance / explained-variance language does not apply the way it does for MOFA: DIABLO's components are
    chosen to discriminate, not to summarise variance.
  - The example data are SYNTHETIC, with a planted structure, so the method can be checked (below). Do not read
    its biology.

Example data and how it was checked
---------------------------------------
Data3/example_mofa_transcriptomics.csv, example_mofa_proteomics.csv and example_mofa_metabolomics.csv (shared with
MOFA — see README-MultiOmics-MOFA.txt for how they were generated): 36 samples, 18 Control / 18 Treated, with a
treatment signal (Reactome's citric acid cycle and glycolysis genes, their proteins, and the matching metabolites)
planted to correlate with the groups, plus a cholesterol and a lipid signal that are real but NOT related to the
groups, plus ~600 noise genes and missing protein values. On these data DIABLO's cross-validated error is close to
0 (its permutation test gives p ≈ 0.05, the smallest a 20-permutation test can report), and it selects the planted
treatment genes (HKDC1, NUP205, ISCA2, SDHA, PFKP, ...) and metabolites (2-Oxoglutarate, Fumarate, Oxaloacetate,
cis-Aconitate, Pyruvate — the TCA-cycle intermediates), not the cholesterol or lipid features. The pathway test of
the selected transcripts finds Glycolysis and Glucose metabolism significant (FDR < 0.05).

Method references
-------------------
  Singh A, Shannon CP, Gautier B, et al. (2019) DIABLO: an integrative approach for identifying key molecular
      drivers from multi-omics assays. Bioinformatics 35:3055-3062.
  Rohart F, Gautier B, Singh A, Lê Cao KA (2017) mixOmics: An R package for 'omics feature selection and
      multiple data integration. PLoS Comput Biol 13:e1005752.
  Tenenhaus A, Philippe C, Guillemot V, et al. (2014) Variable selection for generalized canonical correlation
      analysis. Biostatistics 15:569-583.
  Lê Cao KA, Rossouw D, Robert-Granié C, Besse P (2008) A sparse PLS for variable selection when integrating
      omics data. Stat Appl Genet Mol Biol 7:Article 35.
  Fisher RA (1922) On the interpretation of chi-square from contingency tables, and the calculation of P.
      J R Stat Soc 85:87-94.
  Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. J R Stat Soc B 57:289-300.
  Reactome pathway sets: Milacic M, et al. (2024) The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res
      52:D672-D678.
