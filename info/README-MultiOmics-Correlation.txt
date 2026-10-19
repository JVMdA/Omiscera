Multi-omics — Correlation-based Integration
===========================================

What this method is for
-----------------------
Joint Pathway Enrichment and Cross-omics Pathway Mapping work on lists and
differential-expression tables. This method needs something they cannot use: a
table of VALUES per SAMPLE for each omics, measured on the same samples. It asks
whether the omics layers move together across those samples:

  Same-gene concordance   For every gene present in both the transcriptomics and
                          the proteomics table: how well does its transcript
                          level track its protein level across the samples?
  Cross-omics pairs       Which feature of one omics is correlated with which
                          feature of another (transcript - metabolite,
                          protein - metabolite, transcript - protein)?

IT IS DESCRIPTIVE. A correlation across samples is an association, not a
mechanism (see "What the numbers do not tell you").

How to use it
-------------
1. Open the Multi-omics box "Correlation-based Integration". Choose two or three
   of Transcriptomics, Proteomics, Metabolomics. Genomics cannot be chosen: a
   gene list has no values per sample.
2. Pick the organism (needed to map transcript and protein identifiers to gene
   symbols so a transcript and its protein can be paired).
3. Upload one table per omics (formats below). Each upload tells you what was
   understood: the layout it detected and why, the first sample ids, and any
   column set aside. If the layout is wrong, change "Orientation" - the same file
   is read again.
4. "Resolve identifiers and match samples". The alignment card lists the samples
   every file has, and, per omics, the ones left out.
5. Setup, then "Run correlation".

Input formats
-------------
One table per omics. A row is a feature and a column a sample, or the other way
round - Omiscera detects which:
  - if a text column is named like a sample grouping (Group, Condition, Treatment,
    Phenotype, Cohort, Genotype, Label - MetaboAnalyst's name - or "Sample type" /
    "Sample class") and some of its values repeat, the samples are in the rows and
    that column is set aside (never treated as a feature);
  - otherwise the axis with FEWER labels is taken as the samples (a study has
    fewer samples than features).
Other text columns are set aside too but do not decide anything: a column that
annotates FEATURES (gene_biotype, Description, Protein names, a compound's Class)
says nothing about where the samples are. A plain "Class" or "Type" column is
ambiguous (sample class or compound class?) so it is not used as evidence - if
the guess is wrong, set the orientation yourself.
The first column (or row) holds the ids. Transcriptomics: gene symbol, Ensembl or
Entrez id. Proteomics: UniProt accession (protein groups like "P1;P2" and
isoforms are handled). Metabolomics: the compound name as you wrote it - nothing
here needs a database id, so compounds are correlated by name.
Values should be ALREADY NORMALISED (log2 CPM / TPM / VST for RNA; log-scaled
intensities for protein). Raw counts or raw intensities: tick "Apply log2(x + 1)"
(refused if any value is negative). Spearman, the default, does not care about the
scale; Pearson does.
Repeated feature ids are kept as separate features ("X", "X (2)"); both are still
mapped to their gene through the original id "X". Repeated SAMPLE ids are an error.
Missing values are allowed.

Sample matching
---------------
Sample ids are matched exactly after trimming spaces; an option ignores upper /
lower case. Only samples that EVERY file has are used, in the first file's order;
the ones left out are listed per omics, never silently dropped. Fewer than 5
shared samples is refused (the message says whether the OTHER orientation of a
file would have shared more); fewer than 10 gives a warning, because the p-values
below rest on an approximation that is unreliable at that size.

What is computed
----------------
Features with more than 20% missing values, or with no variation, cannot be
correlated: they are dropped and counted (Overview, "Features used").
  Correlation   Spearman (default): for each pair, keep the samples BOTH
                features have, replace each feature's values by their ranks
                among those samples (ties get the average rank), and take
                Pearson's r of the ranks - the standard pairwise-complete
                Spearman (what R's cor(use = "pairwise.complete.obs", method =
                "spearman") and pandas give). Pearson: r on the values over the
                shared samples. Each pair carries its own n (at least "minimum
                samples per pair", default 5).
  p-value       t = r * sqrt((n - 2) / (1 - r^2)), two-sided, n - 2 degrees of
                freedom: the usual test of "no correlation". It is an
                approximation; with Spearman and few samples it is rough.
  FDR           Benjamini-Hochberg. The three omics pairs are separate families
                (each is corrected over its own pairs); the same-gene test is its
                own family, corrected over the genes tested. The number of pairs
                tested is shown for each family.
  Reported      A cross-omics pair is REPORTED when its FDR <= the FDR cutoff
                (default 0.05) AND |rho| >= the |rho| cutoff (default 0.6). The
                Overview shows how many pairs passed each, and the expected
                number of false ones among the reported (the sum of their
                q-values).
  Same-gene     Features that map to the same gene within one omics (isoforms,
                protein groups) are averaged first; each gene in both tables is
                one test over ALL its features (no variance filter).
  Cross-omics   To keep the number of pairs manageable, each omics contributes
                its most variable features (default 2000, 50 to 5000). Every
                feature of one omics is tested against every feature of the
                other: nothing else is pre-selected.
A pair list is capped at 20,000 rows (a note says so); the strongest by FDR are
kept.

What you get
------------
  Overview               Pairs tested, significant and reported per omics pair;
                         features used; the limitations below; the sample matching.
  Same-gene concordance  Counts, a histogram of rho across the genes, and a table
                         (gene, rho, p-value, FDR, n, the ids that were averaged).
                         Click a row for that gene's values per sample.
  Cross-omics pairs      One table per omics pair. Click a row to see the two
                         features' values per sample (one point per sample, named).
  Heatmap                The features with the most reported partners on each
                         side, clustered (average linkage) so similar profiles sit
                         together. Colour = rho (choose the gradient: red-blue, orange-
                         purple, green-purple, teal-brown, pink-green or blue-
                         yellow-red, all colour-blind-friendly, with a swap
                         option; the choice also colours the histogram and the
                         network lines); a dot marks a REPORTED pair. Drawn at a fixed cell size inside a frame
                         you scroll.
  Network                The strongest reported pairs of every omics pair, one
                         column per omics, rows ordered to cut crossing lines.
                         Line colour = sign, thickness = |rho|.
Every chart has the zoom / pan / autoscale / reset / download tools row (the heatmap and the
network repeat them as buttons above the figure, with Download figure); every table downloads as CSV.

What the numbers do not tell you
--------------------------------
  - GROUPS. If your samples fall into groups (control and treated) and two
    features both differ between the groups, they will correlate even if they are
    unrelated inside each group: a strong correlation can be the treatment effect
    seen twice. A within-group version is a possible later addition.
  - SMALL n. With 12 samples a Spearman rho near 0.6 is common by chance for any
    single pair; only a very strong, clean relation survives the correction over
    thousands of pairs, and the p-value approximation is rough. Treat the
    reported pairs as leads to inspect (click a row), not findings.
  - mRNA - protein. In published studies the typical gene's transcript-protein
    correlation is moderate (Vogel & Marcotte 2012; Liu et al. 2016) because
    protein levels are also set by translation and degradation. A median well
    above zero is expected; low concordance is not a measurement failure.
  - NORMALISATION. Values are used as given. Different scales or batch effects in
    one omics become correlations.
  - MECHANISM. Correlation is not causation, and an RNA - metabolite correlation
    does not say the enzyme makes that metabolite. Use Cross-omics Pathway
    Mapping to look at the pathway context.
  - MISSING VALUES. A pair with missing values is computed on fewer samples (its
    n is shown), so it needs a stronger rho to reach the same p-value. Nothing is
    imputed.

Running time: with no missing values a block of 2,000 x 2,000 features takes a
few seconds. Ranking every pair over its own samples costs more when values are
missing: about 1-2 s at 12 samples, ~15 s at 50 and ~1 min at 100 samples for
2,000 x 2,000 features with 10% missing on both sides (the progress bar shows
it). Lower "features per omics" to go faster.

The example data
----------------
Data3/example_corr_transcriptomics.csv, example_corr_proteomics.csv and
example_corr_metabolomics.csv are three tables for the SAME 12 samples
(Control_1-6, Treated_1-6): genes and proteins in rows with the samples in the
columns, and the metabolite table with samples in the rows and a Group column, so
the orientation detection is exercised by the example itself.
THEY ARE SYNTHETIC, on purpose: no public dataset with matched RNA, protein and
metabolite levels can be bundled, so the structure is planted and known (made by
tests/generate_correlation_examples.py, seeded): 15 TCA-cycle / glycolysis genes
whose transcript and protein share a factor (concordant), 4 whose protein is its
opposite (discordant), CS / IDH3A / MDH2 driven by citrate / isocitrate / malate,
about 270 background genes with independent noise, and about 1% of the protein
cells missing. Citrate, isocitrate and malate also differ between Control and
Treated in the metabolite table, so part of those correlations IS the group
difference - the caveat above. Not every planted gene survives the correction at
12 samples; that is what the method honestly gives, and the cutoffs are not
loosened to show more.

How it was checked
------------------
tests/test_multiomics_correlation.py recomputes every statistic with independent
pure-python code (ranks, Pearson and Spearman with ties and missing values, a
hand-written incomplete-beta p-value, Benjamini-Hochberg) and checks that on
independent random data about 5% of p-values fall below 0.05, uniformly, and that
nothing survives the correction. tests/live_e2e_multiomics_correlation.py reruns
the example through the real server and rebuilds, from the raw CSV files, every
gene's concordance and the COMPLETE transcript - metabolite block (12,558 pairs):
the reported set must be identical; the protein - metabolite block, whose proteins
have missing values, is rebuilt the same way with every pair ranked over its own
samples. The browser tests draw every chart at its real
size and measure overlaps.

References
----------
Spearman C (1904) The proof and measurement of association between two things.
  Am J Psychol 15:72-101.
Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical
  and powerful approach to multiple testing. J R Stat Soc B 57:289-300.
Vogel C, Marcotte EM (2012) Insights into the regulation of protein abundance
  from proteomic and transcriptomic analyses. Nat Rev Genet 13:227-232.
Liu Y, Beyer A, Aebersold R (2016) On the dependency of cellular protein levels
  on mRNA abundance. Cell 165:535-550.
Rohart F, Gautier B, Singh A, Le Cao KA (2017) mixOmics: an R package for 'omics
  feature selection and multiple data integration. PLoS Comput Biol 13:e1005752
  (the standard pairwise-correlation toolkit this can be compared with).
