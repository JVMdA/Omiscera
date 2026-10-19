Multi-omics — Joint Pathway Enrichment
=======================================

What this module is for
-----------------------
You have more than one layer of data on the same biology - say a list of genes
and a list of metabolites, or an RNA-seq table and a proteomics table - and you
want to know whether the layers point at the SAME pathways.

Multi-omics opens a screen with one box per planned method; the unbuilt ones are
greyed out. It will grow one analysis at a time. Today it has two (this file
describes the first; the second has its own file, README-MultiOmics-Mapping.txt):

  Joint Pathway Enrichment
      Input:  2 to 4 omics - Genomics, Transcriptomics, Proteomics,
              Metabolomics - one uploaded file per omics.
      Asks:   for each KEGG / Reactome pathway, do the omics TOGETHER point at
              it? One joint p-value per pathway, with every omics' own evidence
              shown beside it.
      Databases: KEGG and Reactome.

  Cross-omics Pathway Mapping
      Inside ONE pathway, which members did each omics hit, and how do they
      connect (enzyme -> reaction -> metabolite)? Descriptive, no p-value. See
      README-MultiOmics-Mapping.txt.

Planned next, not built yet: network-based integration (STRING), correlation-
based integration, and import of MOFA / DIABLO results.

Nothing is read from the other modules
--------------------------------------
Every omics is uploaded and resolved INSIDE Multi-omics, with that omics' own
resolver (UniProt for genes, mygene.info for transcripts, UniProt ID mapping for
proteins, PubChem/KEGG plus the local cross-reference for compounds). So an
identifier means the same thing here as in that omics' own tile - and re-running
another tile can never silently change what is being combined here. The cost is
that you upload again; the benefit is that the result is always reproducible
from the files you can see.

What to upload
--------------
  Genomics         A gene list: one gene per row, column "Gene Symbol" or
                   "Ensembl Gene ID". The genes you upload ARE the test set.
                   A gene whose identifier is ambiguous or not found is left
                   out and counted - there is no candidate-picking step here and
                   Omiscera never guesses an identifier.
  Transcriptomics  A differential-expression table: an identifier column, log2
                   fold change and p-value (padj optional - computed if
                   absent). A plain gene list is rejected: Omiscera needs the
                   whole table to know what was TESTED.
  Proteomics       A protein-level differential-abundance table: accession,
                   log2 fold change, p-value (padj optional). Same rule.
  Metabolomics     A compound list: one compound per row in a column "Name"
                   (optional InChIKey / PubChem / KEGG / HMDB / ChEBI columns
                   make matching more reliable). No statistics columns.

For the two table inputs Omiscera picks the significant features itself, with
the padj and |log2 fold change| cutoffs you set (defaults 0.05 and 1.0), up- and
down-regulated POOLED into one list: the joint question is "is this pathway
perturbed", and a compound list has no direction. The reference set (the
"universe") is every gene / protein in your table - the field's standard for
expression data - not the whole database. For Genomics and Metabolomics lists
the reference set is the database's own annotated set, exactly as in those
omics' own over-representation tiles.

How the omics are combined
--------------------------
Step 1, per omics. Each omics is tested against every pathway on its own with a
one-sided Fisher exact test - the same test, and the same reference-set rules,
that omics' own over-representation tile uses. The per-omics p-value you see in a
joint result is the SAME number that tile reports for the same input; this is
checked automatically against the Genomics and Metabolomics tiles (and, offline,
against the Transcriptomics, Proteomics and Metabolomics engines).

Step 2, combine. The per-omics p-values for a pathway are combined into one:

  Brown's method (default) - Fisher's test, corrected for correlated omics
      Fisher's statistic X = -2 * sum(ln p_i) assumes the omics are independent.
      They often are not: RNA and protein from the same samples, or a gene list
      derived from the same experiment, report largely the SAME signal, and
      Fisher counts it twice. Brown (1975) keeps X but matches its mean and
      variance to a scaled chi-squared: with k omics and pairwise correlations
      rho_ij, Var[X] = 4k + 8 * sum(rho_ij), c = Var / 4k, f = 8k^2 / Var, and
      p = P(chi-squared_f >= X / c). The correlations are estimated from this
      run's own p-values (the "empirical" version, Poole et al. 2016): for each
      pair of omics, the Pearson correlation of -2 ln p across every pathway that
      both could test - the whole jointly tested population, not only the
      pathways with a hit. That works even for omics that share no features
      (genes vs metabolites), which a feature-level correlation could not.
      Three properties are tested: with every correlation 0 it IS Fisher; an
      omics combined with an exact copy of itself returns its own p-value (Fisher
      would claim a far smaller one - the double counting); and wherever Fisher's
      p is below about 0.28 (every result that could be significant) Brown's p is
      never smaller than Fisher's. Above ~0.28 the two tails cross and Brown can
      be slightly smaller - a difference between two clearly non-significant
      p-values. A negative estimated correlation is set to 0 (a negative value
      would make Brown claim MORE evidence than independent omics). With fewer
      than 20 jointly tested pathways the correlation is too noisy to estimate:
      the omics are then treated as independent (plain Fisher) and the result SAYS
      so. The results page lists every pair's correlation, and how many
      independent omics your layers amount to.
  Fisher's combined probability test
      X2 = -2 * sum(ln p_i), compared with chi-squared on 2k degrees of freedom
      (k = number of omics). Evidence ACCUMULATES: several omics each moderately
      enriched can add up to a strong joint result - but it assumes independence,
      so correlated omics give a joint p-value that is too small. Kept for
      comparison; a Fisher run still REPORTS the estimated correlation, marked as
      information only. (R.A. Fisher, 1932, "Statistical Methods for Research
      Workers".)
  Tippett's minimum-p test
      p = 1 - (1 - min p)^k. Asks whether AT LEAST ONE omics is strongly
      enriched in the pathway, ignoring agreement between them. (L.H.C. Tippett,
      1931, "The Methods of Statistics".)

Merging per-dataset p-values, with Brown's correction for correlated datasets, is
also how the published tool ActivePathways combines omics (Paczkowska et al.
2020, Nature Communications 11:735). This is NOT the "pool every hit into one
bigger table and run one test" approach: combining each omics' own p-value keeps every omics'
evidence separate and visible, which a pooled count would hide. The code is
Omiscera's own; the combination arithmetic itself is scipy's
scipy.stats.combine_pvalues, an already-installed, correct implementation, so it
is not re-written by hand.

Stouffer's Z method is deliberately NOT offered. It treats a p-value near 1 as
strong evidence AGAINST a pathway (z = -infinity at p = 1), but a one-sided
over-representation p-value of 1.0 only means "no enrichment seen" - and an omics
with no hit in a pathway legitimately has p = 1.0. Measured: one omics at
p = 1e-10 plus one at p = 1.0 combines to p = 0.90 under Stouffer (a real hit
wiped out) but 2.4e-9 under Fisher. Brown, Fisher and Tippett are unaffected by p = 1.

The rules that decide what is combined (each one tested)
--------------------------------------------------------
  1. An omics with no hit in a pathway enters as p = 1.0. The single-omics tiles
     drop zero-hit pathways from their own correction (right for a standalone
     table); the combination must not, because "this omics saw nothing here" is
     real evidence when another omics did.
  2. A pathway enters the combined family only if SOME omics has a hit. A pathway
     nobody hit says nothing, and counting it would only inflate the correction.
  3. A pathway is combined only if it passes the size gate in EVERY usable
     omics. One that passes in only some is listed in a separate "Partial
     coverage" table, with each omics' own result where it exists - never padded
     with a made-up p-value.
  4. An omics with none of its hits annotated in a database cannot be tested
     there (its p-values would be undefined, not 1.0): it is left out of that
     database's combination, with a note. Fewer than two usable omics gives an
     empty result and says why, never a one-omics "combination".
  5. Benjamini-Hochberg correction is applied ONCE, to the combined p-values
     only. The per-omics p-values beside them are informational. Each database is
     corrected on its own; there is no correction across databases.
  6. KEGG's organism is resolved ONCE per run (the same shared resolver the
     Transcriptomics and Proteomics tiles use) and the same code is handed to
     every omics' KEGG builder. Resolving it separately per omics once let two
     code paths disagree on the organism and silently test the wrong maps. If no
     KEGG organism can be resolved (no genome, or several strains and none
     confirmed) the KEGG result is EMPTY and says why - never a fallback to
     generic maps. Reactome is aligned through the same organism name on every
     omics.

Why KEGG and Reactome only
--------------------------
A joint test needs a database that annotates BOTH genes/proteins AND compounds
under the same pathway identifiers. KEGG and Reactome do, and both are already
built for every omics here. WikiPathways exists in Omiscera only as a gene-only
collection (no compound side yet), and MetaCyc's bulk files are likely
licence-restricted - both are candidates for later, not v1.

Reading the results
-------------------
  Combined table   One row per jointly tested pathway: each omics' hits (k), its
                   annotated members (K), its own p-value and the overlapping
                   members, then the combined p-value, combined FDR, the
                   statistic, the method, and how many omics were combined.
  Partial coverage Pathways tested in only some omics (see rule 3). No joint
                   p-value.
  Joint dot plot   One row per pathway, most significant first. Position =
                   -log10 of the combined p-value; colour = -log10 of the combined
                   FDR; circle size = hits summed over the omics.
  Evidence by omics
                   Which omics drove each pathway: one coloured circle per omics
                   at -log10 of its own p-value, and a black diamond at the
                   combined one. A diamond far to the right of every circle marks
                   a pathway no single omics shows on its own - that is what
                   combining adds. A circle at 0 means that omics saw nothing
                   there. Each pathway row is split into one thin lane per
                   omics (and one for the diamond), so markers never hide each
                   other even when their p-values are identical.

Limits to keep in mind
----------------------
  - Brown's correlation is ESTIMATED from this run's own p-values, so it carries
    sampling noise and also picks up real shared biology (pathways that are truly
    active in several omics look "correlated"). That makes the correction
    somewhat conservative. Below 20 jointly tested pathways it is not estimated
    and the omics are treated as independent (plain Fisher), with a note.
    Plain Fisher and Tippett do not correct for correlation (Tippett is already
    conservative under positive dependence). Either way, read the combined
    p-value and FDR as a RANKING of pathways, not as an exact error rate.
    (Every result says so.)
  - Transcriptomics here uses the plain hypergeometric test. That module's own
    optional gene-length-bias correction is not applied, so a joint p-value can
    differ slightly from a Transcriptomics run that used it.
  - One organism for the whole run. Omics measured in different species cannot
    be combined here.
  - Only the most recent uploads and resolutions of the current session are kept.
    Re-uploading an omics clears the resolution and any result built on it.
  - KEGG covers only organisms with a sequenced, annotated genome; Reactome's
    curated (non-electronic) coverage is essentially human-only.
  - Genes that are ambiguous or not found in Genomics are left out (counted on
    the resolve screen and in the results).

Downloads
---------
Each database's combined table and partial-coverage table download as CSV with
the full, unclamped member lists (the on-screen cells are clamped to three lines
with the full text in the tooltip, so one long list cannot blow up a row).

Example files
-------------
Data3/example_multiomics_genes.csv        29 human TCA-cycle / glycolysis genes
Data3/example_metabolomics_list.csv       TCA-cycle / glycolysis compounds + a few lipids
Data3/example_transcriptomics_de_full.csv a full transcriptome DE table
Data3/example_proteomics_de_full.csv      a full proteome DE table
The gene and compound examples plant the same biology on purpose: the joint
result must put KEGG's TCA cycle and glycolysis at the top, with a combined
p-value stronger than either omics alone (checked live; TCA cycle: 7e-50 jointly
against 6e-41 for the better single omics).

Citations
---------
  Fisher RA (1932). Statistical Methods for Research Workers, 4th ed. Oliver and
      Boyd, Edinburgh.
  Brown MB (1975). A method for combining non-independent, one-sided tests of
      significance. Biometrics 31:987-992.
  Poole W, et al. (2016). Combining dependent P-values with an empirical
      adaptation of Brown's method. Bioinformatics 32:i430-i436.
  Tippett LHC (1931). The Methods of Statistics. Williams and Norgate, London.
  Paczkowska M, et al. (2020). Integrative pathway enrichment analysis of
      multivariate omics data. Nature Communications 11:735.
  Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. J R Stat
      Soc B 57:289-300.
