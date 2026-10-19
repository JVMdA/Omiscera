Multi-omics — MOFA / Factor-based Integration
=============================================

What this method is for
-----------------------
The other Multi-omics methods start from lists, differential-expression tables
or single pairs of features. MOFA (Multi-Omics Factor Analysis) asks a different
question of a table of VALUES per SAMPLE for each omics, measured on the same
samples: which few hidden FACTORS explain the differences between the samples,
and how do they spread over the omics?

For every factor MOFA gives you:
  Variance explained   How much of each omics' variation the factor accounts
                       for. A factor that explains a lot in several omics is a
                       SHARED pattern (the omics move together); one that
                       explains variance in a single omics is specific to it.
  Factor values        One number per sample: where the sample sits on the
                       factor. This is what you plot, colour by group, or
                       correlate with a phenotype.
  Weights              One number per feature: how strongly the feature belongs
                       to the factor. The top-weighted features (and the
                       pathways they point to) tell you what the factor is about.

IT IS DESCRIPTIVE. A factor is a pattern of co-variation, not a pathway and not
a cause (see "What the numbers do not tell you").

How it is run
-------------
Omiscera runs the AUTHORS' OWN implementation, the Python package mofapy2
(version 0.7.5, LGPL-3), as a library, with a Gaussian likelihood in every omics.
It does not re-implement the model. What Omiscera adds is around it: the data
preparation, a fixed order and sign for the factors, a stability check, the
group tests and the pathway test of the weights - all written out below and
tested against independent calculations (tests/test_multiomics_mofa.py).

The project's Python environment matters: mofapy2 needs a working numpy / scipy.
Use the virtual environment described in requirements.txt (the Spyder runtime on
the development machine could not run any matrix product at all).

How to use it
-------------
1. Open the Multi-omics box "MOFA / Factor-based Integration". It opens the same
   upload screen as Correlation-based Integration: choose two or three of
   Transcriptomics, Proteomics, Metabolomics (Genomics cannot take part - a gene
   list has no values per sample), pick the organism, upload one table per omics,
   and "Resolve identifiers and match samples". Formats and the way the samples
   are matched are exactly those of Correlation-based Integration (see
   README-MultiOmics-Correlation.txt). The "Example file" links give a synthetic
   36-sample example for each omics (the Correlation examples have 12 samples,
   too few for MOFA).
2. "Continue to MOFA setup": number of factors, convergence, features kept,
   scaling, stability refits, and the sample groups (below). "Run MOFA".
3. Results panes: Overview, Variance explained, Factors, Weights, Pathway
   enrichment (run on request, per omics and factor), Data used, Parameters used.
   Every figure has the zoom / pan / in / out / autoscale / reset buttons, a Size
   menu and "Download figure" (SVG / PNG / TIFF); every table can be downloaded
   as CSV.

Sample size
-----------
MOFA needs at least 6 samples (fewer is refused). Below about 15 the factors, and
above all the weights, are unreliable; the setup and the results say so. The
authors' own guidance is to use tens of samples or more.

Data preparation
----------------
For every omics, in this order:
  - optional log2(x + 1) (for raw counts / intensities; refused for negative
    values);
  - features missing in more than 50% of the samples are dropped; constant
    features are dropped. Missing values that remain are LEFT missing - MOFA
    handles them, nothing is imputed by Omiscera;
  - the N most variable features are kept (default 2000; MOFA's authors advise
    using highly variable features);
  - features are centred, and, with "Give every omics the same total variance"
    (default on, MOFA's scale_views), each omics is scaled so the omics with the
    most features does not dominate.

The model
---------
MOFA is a Bayesian factor model: each omics' data matrix (samples x features) is
written as (factor values) x (weights) + noise, with sparsity priors so that a
factor is switched off in the omics that do not need it. You offer a number of
factors; MOFA drops the ones the data do not support, and Omiscera additionally
drops any factor explaining less than a minimum (default 1%) in EVERY omics.

OFFER MORE FACTORS THAN YOU EXPECT. This is the authors' advice, and we measured
why on the bundled example (three planted factors): offering 3, 4, 5 or 6 factors
lost one of the real factors every time (a different one each time); offering 8
or more recovered all three. The default is 10.

Order and sign. MOFA's factors come out in an arbitrary order and sign. Omiscera
orders them by the total variance they explain (largest first) and flips each so
that its largest weight (over all omics) is positive. A re-run of the same data
and settings therefore gives the same table. The sign has no biological meaning:
a factor can be read upside down.

Variance explained uses mofapy2's own definition: R2 = 1 - (residual sum of
squares) / (total sum of squares), per omics and factor, on the data the model saw
(missing cells enter as the model's own imputed values with the prediction set to
zero there). The test suite recomputes it from the factor values and weights with
an independent implementation and gets the same numbers to 1e-8.

Is the result stable?
---------------------
MOFA starts from a random point and can find a different local optimum. Its
authors advise refitting several times and trusting the factors that come back.
Omiscera refits with other seeds (default 3 refits, seed 42 for the main fit);
for each factor it reports the median, over the refits, of the best absolute
correlation with any factor of that refit. Below the cut-off (default 0.7) the
factor is flagged UNSTABLE: shown, never hidden, and the warning says to read it
as noise, not as a finding.

Sample groups (interpretation only)
-----------------------------------
The model is fitted on ALL samples together; groups play no part in the fit. They
colour the plots and let each factor be tested for a difference between groups.
Where do the labels come from?
  - a label column of a table (the Group column of a metabolite table with
    samples in the rows is offered automatically), or
  - a small sample sheet you upload: the first column is the sample id, every
    other column is a label (Group, Batch, Sex ...); a column is offered when it
    has a few repeated values.
Test per factor: Mann-Whitney U for two groups (effect = rank-biserial correlation,
-1 to +1), Kruskal-Wallis for more (effect = eta squared from H, (H - k + 1) / (n - k), 0 to 1;
Tomczak & Tomczak 2014); groups with
fewer than 2 samples are left out; Benjamini-Hochberg across the factors. The
p-values equal scipy's, the FDR equals the one used by every other module.
Remember: if the groups differ in the way you designed them (control vs treated),
the factor that follows the groups is the treatment effect seen across the omics.
That is expected, and it is not a discovery.

What the weights tell you
-------------------------
A large positive weight: the feature is high in the samples with a high factor
value; a large negative weight: it is low there. Only the sign RELATIVE to the other
features of the same factor means something. The Weights pane shows the 25 largest
|weights| of a chosen factor and omics as a bar chart (two features that carry one
gene symbol get their feature id in the label so both stay visible) and a table; the
CSV has every feature of the omics for every factor.

Pathway enrichment of a factor's weights
----------------------------------------
Which Reactome pathways does a factor point to? For one omics and one factor, the
weights of the pathway's members are compared with the weights of all the other
features: a Mann-Whitney U test (continuity- and tie-corrected, equal to scipy's
asymptotic test), the idea of parametric / rank gene-set analysis of loadings
(PCGSE, Frost et al. 2015) which MOFA's own enrichment function follows.
  Signed          higher or lower weights? (AUC above 0.5 = the members have
                  higher weights; below = lower). Two-sided.
  Absolute        larger |weights|, whatever the sign? One-sided.
Benjamini-Hochberg across the pathways tested. Pathways are Reactome sets built
with the same builders as Joint Pathway Enrichment: transcripts by Entrez Gene id,
proteins by UniProt accession, metabolites by matched compound (the compound-name
matching and the pathway sets are built once per result, on the first request for that
omics: seconds for transcripts and proteins, several minutes for metabolites; later
factors and tests reuse them). Features that match
no pathway are left out of the test; several features of one gene count once (the
one with the largest |weight| represents it). Pathways with fewer than a minimum
(default 10) or more than a maximum (default 500) matched members are not tested.
Run per omics and factor, on request, so the pathway sets are built only when you
ask.

What the numbers do not tell you
--------------------------------
  - A factor is not a pathway and not a cause. Read its biology from the top
    weights and the pathway test, and take it as a lead to test.
  - Variance explained is a share of the variation of the features that entered
    the model, not of "the biology".
  - With few samples the weights are noisy; with many features and few samples,
    some factors will be driven by one or two samples - look at the factor values.
  - The sign and the order of a factor are conventions.
  - A factor that follows your sample groups is the group effect seen across the
    omics.
  - The example data are SYNTHETIC: the factors were planted so that the method
    can be checked (below). Do not read their biology.

Example data and how it was checked
-----------------------------------
Data3/example_mofa_transcriptomics.csv, example_mofa_proteomics.csv and
example_mofa_metabolomics.csv: 36 samples (18 Control, 18 Treated), made by
tests/generate_mofa_examples.py (seeded; needs the network once for gene symbols).
Three factors are PLANTED: a treatment factor shared by all three omics (the genes
of Reactome's citric acid cycle and glycolysis, their proteins and the TCA /
glycolysis metabolites); a cholesterol-synthesis factor in transcripts and
proteins only (unrelated to the groups); a lipid factor in metabolites only
(unrelated to the groups); 600 noise genes; about 4% missing protein values. The
true factor values are in tests/fixtures/mofa_example_truth.csv. On these data
MOFA recovers the factors with |r| of 0.998 (treatment), 0.990 (cholesterol) and
0.987 (lipids); the treatment factor is the only one that differs between the
groups; the fourth factor it keeps (1% of the variance) is flagged unstable; and
the pathway test of the treatment factor puts Glycolysis and Glucose metabolism
first.

Method references
-----------------
  Argelaguet R, Velten B, Arnol D, et al. (2018) Multi-Omics Factor Analysis - a
      framework for unsupervised integration of multi-omics data sets. Mol Syst
      Biol 14:e8124.
  Argelaguet R, Arnol D, Bredikhin D, et al. (2020) MOFA+: a statistical framework
      for comprehensive integration of multi-modal single-cell data. Genome Biol
      21:111.
  Frost HR, Li Z, Moore JH (2015) Principal component gene set enrichment (PCGSE).
      BioData Min 8:25.
  Mann HB, Whitney DR (1947) On a test of whether one of two random variables is
      stochastically larger than the other. Ann Math Stat 18:50-60.
  Kruskal WH, Wallis WA (1952) Use of ranks in one-criterion variance analysis.
      J Am Stat Assoc 47:583-621.
  Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. J R Stat
      Soc B 57:289-300.
  Tomczak M, Tomczak E (2014) The need to report effect size estimates revisited. Trends
      Sport Sci 21:19-25.
  Reactome pathway sets: Milacic M, et al. (2024) The Reactome Pathway Knowledgebase
      2024. Nucleic Acids Res 52:D672-D678.
