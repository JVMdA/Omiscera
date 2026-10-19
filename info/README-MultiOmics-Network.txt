Multi-omics — Network-based Integration
=======================================

What this method is for
-----------------------
Joint Pathway Enrichment asks "which pathways do all my omics point at?". Cross-omics
Pathway Mapping shows where, inside one pathway, each omics lands. This method asks
two questions that need a NETWORK:

  1. THE STRING NETWORK (genes / proteins). Which proteins sit at the centre of the
     COMBINED response of all your omics - including proteins nobody measured
     ("recruited") - and which connected modules do the hits form? Each module shows
     which omics fed it.
  2. THE REACTOME PATHWAY NETWORK (pathways). Which of the pathways your omics hit
     overlap in the members that were hit, so that they form groups of related
     pathways?

IT IS DESCRIPTIVE. A protein's score is a position in a network built from literature
and databases, not a measurement; a module is a group of well-connected proteins that
carry signal, not a tested hypothesis.

How to use it
-------------
1. Open the Multi-omics box "Network-based Integration". Choose two to four omics,
   upload one file each, pick the organism, resolve the identifiers (exactly as for
   Joint Pathway Enrichment; if you already resolved your files in this session the box
   takes you straight on - nothing is uploaded twice).
2. Setup: the padj / |log2 fold-change| cutoffs for the differential-expression tables,
   the STRING confidence cutoff, the number of permutations, the significance threshold,
   the module settings, the metabolite reaction limit and the size of the pathway network.
3. Run. The first run downloads the STRING network and Reactome's reaction files (kept
   afterwards); the permutation test is the slow step (about 0.3 s per permutation on the
   human network).

STRING needs a network for your organism (about 12,500 species have one); if STRING has
none, the box says so by name and the other Multi-omics methods still work.

The STRING network
------------------
Seeds. Every omics supplies seeds, proteins placed on the STRING network:
  Genomics         the genes of your list.
  Transcriptomics  the genes passing the padj and |log2FC| cutoffs, up- and down-
                   regulated pooled - the SAME selection the Transcriptomics tiles use.
  Proteomics       the proteins passing the cutoffs, by the Proteomics tiles' own rule.
  Metabolomics     STRING is a protein network, so a hit metabolite seeds the proteins
                   that Reactome lists in the reactions it takes part in (below).
Entities that cannot be placed on STRING are counted and shown, never guessed.

Equal weight per omics. Each omics' seeds are scaled to sum to 1, the scaled values are
averaged over the omics that have seeds, and the result is multiplied by the number of
seeds so the mean seed weighs 1.0:
    weight(protein) = (number of seeds) x (1 / K) x sum over omics of [protein seeded by
                      that omics] / (number of seeds of that omics)
So an omics with 3,000 hits does not drown one with 40, and a protein hit by several
omics gets a larger seed. With ONE omics every seed weighs exactly 1.0 - the result is
the same as that omics' own network tile (checked by a test that compares every score).

Propagation. A random walk with restart from the seeds (Koehler et al. 2008) gives each
protein a score; each protein's score is compared with a degree-preserving permutation
null (the number of permutations is a setup option) to give a Z-score, a p-value and a
Benjamini-Hochberg FDR. These are the numerics of the Genomics / Transcriptomics /
Proteomics network tiles, imported unchanged.

What the hit table means.
  Seed        a protein one of your omics supplied. A seed that scores significantly is
              largely EXPECTED: it receives signal directly and is compared with
              permutations in which it is usually not a seed.
  Recruited   a significant protein no omics supplied. These are the discovery, judged
              against a family of non-seed proteins only (a stricter test), and are the
              default view of the hit table.

Modules. The significant proteins are clustered with MCL (van Dongen 2000) and each
module is labelled with its most enriched GO / Reactome / Complex Portal term. The
columns 'Seeds from <omics>' count how many members each omics supplied, and a bar chart
shows it. A module fed by ONE omics is that omics' own story; one fed by several is where
they agree. The label is a description chosen after the modules were found; its FDR is
not a hypothesis test.

Metabolites through their enzymes
---------------------------------
Reactome publishes which proteins and which compounds take part in each reaction. A hit
metabolite is matched to ChEBI (and, when Reactome cites another protonation state of the
same compound, to the same InChIKey skeleton - as in Cross-omics Pathway Mapping), and
seeds the proteins listed in its reactions.
  - HUBS ARE LEFT OUT. ATP, water, NAD+, CoA or H+ sit in hundreds of reactions and would
    seed a large part of the proteome. A compound in more reactions than the limit
    (default 60 - for human data about 2% of compounds are above it) is not used; the
    'Metabolite enzymes' pane lists every compound, its reactions, its proteins and, if
    it was left out, why. The limit is a setup option.
  - FAMILY REACTIONS ARE SKIPPED. Reactome lists every receptor of a family in a ligand-
    binding step ("Liganded Gi-activating GPCRs bind Gi" lists 212 proteins; the E3
    ubiquitin steps list 254), so a metabolite that is a ligand would seed hundreds of
    unrelated proteins from one step (measured: lactate seeded 298 proteins from 20
    reactions, caffeine 217 from 8). A reaction that lists more than 30 proteins (a setup
    option; 959 of the 15,950 human reactions are above it) does not seed; the metabolite
    table counts these per compound.
  - A LISTED PROTEIN IS NOT NECESSARILY THE CATALYST. Reactome's files list every protein
    taking part in a reaction (enzyme, complex member, regulator, substrate). Read a
    metabolite's proteins as "proteins in the reactions of this metabolite", an
    association.

The Reactome pathway network
----------------------------
Joint Pathway Enrichment's own per-omics tests are run on Reactome and combined (Brown's
method). The pathways with the best combined p-values (a setup option, default 40) that
have at least one hit become nodes - after NEAR-DUPLICATES are merged: Reactome pathways
sit in a hierarchy, so a parent and its sub-steps often share almost the same hits and would fill
the picture (measured on the bundled examples: the 40 best pathways gave 274 links; they are only
16 distinct pathways). Walking the pathways in order of combined p-value, one whose hit members
overlap an already drawn pathway by a Jaccard of at least 0.7 (a setup option) is merged into the
better-ranked one - the hover and the CSV list what it stands for - and the walk goes on until
the chosen number of DISTINCT pathways is drawn. Two pathways are joined when the Jaccard similarity of
their HIT members (genes pooled across omics by symbol; metabolites kept apart) reaches
the threshold (default 0.4, close to the 0.375 Merico et al. used) - the idea of the Enrichment Map (Merico et al. 2010). Colour
= a group of connected pathways; the chips under a pathway name the omics that hit it;
size = number of hit members; a black outline = combined FDR <= 0.05. It is drawn at a
fixed size inside a frame you scroll. Above the figure sit the same options every other figure
has: Size (fit to page, 50% to 200%), Pathway names (all, the 10 best, none), Download figure
(PNG / SVG / TIFF), and the Zoom, Pan, Zoom in / out, Autoscale and Reset axes tools - repeated
as buttons because the chart's own toolbar sits at the far corner of a wide figure. The layout is computed in pixels with a label box around every pathway (its
name over up to three lines, the marker and its chips), and the canvas grows until no two
boxes overlap - so no two labels can cover each other, however many pathways there are.

What the numbers do NOT tell you
--------------------------------
  - STRING's combined score includes text-mining and cross-organism transfer, so it is not
    pure experimental evidence and it favours well-studied proteins. A recruited protein
    is often a well-known one for that reason.
  - A module's label is descriptive, not a test. A module is chosen for being connected and
    carrying signal, then labelled.
  - Equal weight is a choice, not a fact: it says "each omics speaks with one voice". An
    omics that is noisy or has very few hits gets the same total weight as a strong one.
  - Metabolite seeds are an association through Reactome, and only compounds Reactome
    knows can take part.
  - The p-values come from a permutation null that is skewed; trust the ranking more than
    any single FDR. With few seeds or a small network nothing may reach significance - the
    result then says so rather than relaxing the threshold.
  - Joining pathways by shared hit members shows which pathways are the same biology seen
    twice (many pathways share genes); it does not say they interact. Merging near-
    duplicates is a choice to make the picture readable: set the merge threshold to 1 to
    merge only pathways with IDENTICAL hit members; the CSV keeps what was merged.

The example data
----------------
The bundled Data3 examples (a gene list, two differential-expression tables and a
metabolite list) are synthetic: they plant TCA-cycle and glycolysis genes and metabolites in
the Genomics and Metabolomics lists, not in the two expression tables.
What the real run gives on them (human STRING at confidence 0.7, 500 permutations, the
default limits): 1,126 seeds - Genomics 29, Transcriptomics 552, Proteomics 115, and 467
proteins from the 40 usable metabolites (4 compounds left out as hubs: 2-oxoglutarate,
succinate, palmitic acid, glycine; 1 not in Reactome). Eight proteins are significant and ALL
eight are seeds: no recruited protein and no module. That is the honest answer for a small,
diffuse seed set and the thresholds are not loosened to show more. The planted TCA proteins
(CS, SDHA, PDHB, MDH2, FH, SDHB, DLD, ACO2) are among the ten highest propagated scores.
In the pathway network the synthetic proteomics example puts about 60 hits in ribosome and
translation pathways, which take the best combined p-values; with the near-duplicates merged
the network shows the distinct pathways around them (the hover says which pathways each one
stands for), and Reactome's citric acid cycle is the third of them. A large group of cell-cycle
and DNA-repair pathways that genuinely share the same proteomics hits stays densely linked;
raise the link threshold in the setup to thin it.
A run takes about 6 minutes at 500 permutations; the FIRST run also looks up each of
Reactome's ~3,200 compounds at ChEBI (a few minutes, remembered afterwards).

How it was checked
------------------
tests/test_multiomics_network.py checks the equal-weight rule against hand computations,
the metabolite - enzyme join against a brute-force join, the pathway network against a
brute-force Jaccard, and that identical seeds give the existing network tile's scores
exactly. tests/live_e2e_multiomics_network.py runs the whole example on the real STRING
network and re-derives from the raw files the seeds of every omics, the metabolite -
enzyme join (Reactome's files streamed again), the combined weights, THE RANDOM WALK ITSELF
(every protein's score on a graph rebuilt from STRING's links file, to 4e-10) and the
pathway edges. The browser test draws every pane at its real size.

References
----------
Szklarczyk D et al. (2023) The STRING database in 2023: protein-protein association
  networks and functional enrichment analyses for any sequenced genome of interest.
  Nucleic Acids Res 51:D638-D646.
Koehler S, Bauer S, Horn D, Robinson PN (2008) Walking the interactome for prioritization
  of candidate disease genes. Am J Hum Genet 82:949-958.
van Dongen S (2000) Graph clustering by flow simulation. PhD thesis, University of Utrecht.
Merico D, Isserlin R, Stueker O, Emili A, Bader GD (2010) Enrichment Map: a network-based
  method for gene-set enrichment visualization and interpretation. PLoS ONE 5:e13984.
Milacic M et al. (2024) The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res
  52:D672-D678.
Brown MB (1975) A method for combining non-independent, one-sided tests of significance.
  Biometrics 31:987-992.
