DIABLO reference fixtures (mixOmics), for test_multiomics_diablo.py
=====================================================================

Omiscera's own Python implementation of DIABLO is checked here against the published algorithm's own R package,
mixOmics, run once on a fixed, small dataset. The shipped fixtures (tests/fixtures/diablo_reference/) are the
result of that run, committed to the repo; running the OFFLINE test suite needs no R at all - only regenerating
the fixtures does.

To regenerate (only needed if Multiomics_diablo.py's algorithm, or the fixture data, changes):

1. Install R (4.3 was used here) and mixOmics:
     install.packages("BiocManager")
     BiocManager::install("mixOmics")

2. From the repo root, rebuild the input tables and the fixed train/test split (pure Python, offline):
     python tests/generate_diablo_r_fixtures.py

3. Run the real mixOmics fit and prediction on those exact tables:
     Rscript tests/r_reference/diablo_reference.R

   This writes tests/fixtures/diablo_reference/out/*.csv: the training variates and loadings, the AVE (variance
   explained) per block and component, the per-block weights, the selected features, and - on the held-out test
   split - the predicted scores and the class call of every one of mixOmics' three distances (max.dist,
   centroids.dist, mahalanobis.dist), both per-omics and combined (majority vote, weighted vote).

4. Re-run tests/test_multiomics_diablo.py (the usual Python interpreter, no R) and confirm it is still green.

Data: the 36-sample Multi-omics MOFA example tables (Data3/example_mofa_*.csv - Control/Treated, missing protein
values), prepared exactly as Multiomics_diablo.prepare_blocks would (top 60 most-variable transcripts/proteins,
all 42 metabolites). design = 0.1 between every pair of omics, 1.0 from every omics to the outcome (mixOmics'
own convention when a 3x3 design is extended to include Y); ncomp = 2; keepX = 10 then 8 features per omics per
component (small numbers, on purpose: enough to be non-trivial, small enough that the whole run - fit, AVE,
selected features, and prediction under all three distances and both vote rules - takes seconds, not minutes).
