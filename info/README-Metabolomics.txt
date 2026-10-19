Metabolomics — five enrichment analyses
========================================

What this module is for
-----------------------
You have measured (or identified) a set of metabolites and you want to know
what they have in common: do they concentrate in particular pathways, or in
particular chemical or lipid classes - and, for a pathway hit, is it a real
biochemical bottleneck or a peripheral one?

Metabolomics is five separate analyses, one tile each. They answer different
questions, so they are deliberately not folded into one screen:

  Over-Representation Analysis
      Input:  a list of compound names.
      Asks:   does your list contain more members of a PATHWAY than chance
              would give?  Databases: KEGG, Reactome, SMPDB, PathBank.

  Chemical Class Enrichment
      Input:  a list of compound names.
      Asks:   is your list enriched in a CHEMICAL class (e.g. "TCA acids",
              "Xanthines")?  Database: RefMet's taxonomy.

  Lipid Class / Structural Enrichment
      Input:  a list of compound names.
      Asks:   is your list enriched in a LIPID class (e.g.
              "Glycerophosphocholines")?  Database: LIPID MAPS.

  Quantitative Enrichment
      Input:  a concentration (or intensity) table with a group label.
      Asks:   do the compounds in a term shift up or down between two groups,
              compared with all the other compounds you measured?  Any of the
              databases above.  Omiscera computes the per-compound statistics
              for you.

  Pathway Topology Analysis
      Input:  EITHER a list of compound names OR a concentration table -
              your choice, made on the tile itself.
      Asks:   beyond whether a KEGG pathway is enriched, WHERE do your
              matched compounds sit in that pathway's own reaction network -
              are they hub/bottleneck compounds, or peripheral ones?
              Database: KEGG only.

Pathway Topology Analysis in detail
------------------------------------
This is the MetaboAnalyst/MetPA-style "pathway impact" score: a topology
measure of each compound's position in the pathway's reaction graph,
combined with the same enrichment p-value the other KEGG-testing tiles
already compute.

  Impact = (sum of your matched compounds' importance) /
           (sum of every compound's importance in that pathway)

  Importance is either:
    Relative betweenness centrality (default) - how often a compound lies
      on the shortest reaction path between two others. A compound that
      sits between many other pairs is more central to the pathway.
    Out-degree - how many compounds it directly feeds via a reaction.

  Both are computed by Omiscera's own implementation of Brandes' algorithm
  (Brandes, 2001, J. Math. Sociol. 25:163-177) - no third-party graph
  library. The impact definition follows Xia & Wishart, 2010, Bioinformatics
  26:2342 (the paper behind MetaboAnalyst's MetPA tool); topology-aware
  pathway analysis as an approach more generally follows Tarca et al., 2009
  (SPIA). This is Omiscera's own code from these published definitions, not
  a copy of any other tool's implementation.

  The ENRICHMENT test (the p-value/FDR columns) is computed by calling the
  exact same code Over-Representation Analysis and Quantitative Enrichment
  use for their own KEGG results - so for the same input and settings, the
  p-values are identical whichever tile you ran. Pathway Topology Analysis
  adds the Impact column on top; it does not change how significance itself
  is computed. This also means it inherits the organism-specificity of KEGG
  described below (Databases and licences -> KEGG), and the same "no genome,
  no fallback" honesty.

  What is NEW here, and does vary by organism: the pathway's own REACTION
  GRAPH (which compound reacts with which). Verified live against KEGG
  (2026-09-29): the human, yeast and E. coli versions of the TCA cycle map
  (hsa00020 / sce00020 / eco00020) share the same 23 compounds but have 22,
  21 and 21 reactions respectively - different edges connecting the same
  nodes. So the reaction graph is fetched and cached per organism, separately
  from the compound-membership data the other KEGG-testing tiles share.

  A pathway KEGG draws with no defined reaction graph (some maps show
  complexes or regulatory relationships rather than a reaction network) has
  its Impact reported as 0 and explicitly flagged in the "Topology note"
  column - never silently indistinguishable from a genuine zero-impact hit.

  The "Pathway impact" plot lists every tested pathway's NAME on the y-axis
  (one row per pathway, same convention as this app's other dot plots),
  Impact on the x-axis, and is autoscaled on first draw so the largest point
  is never cut off at the plot edge. A bigger circle means more of your
  matched compounds were found in that pathway's own graph (this is
  independent of Impact and of significance - a small pathway can score a
  high Impact from very few matches, a large one a low Impact from many).
  Colour = statistical significance (p-value or FDR, your choice). "Terms
  shown" narrows the list to the most significant N if it runs long; it
  shows every tested pathway by default, since Impact is descriptive of
  topology, not a threshold result.

  The "Pathway graph" figure on the results screen draws the selected
  pathway at KEGG's own map coordinates, with your matched compounds
  highlighted in red. Nodes are labelled with the compound's real name
  (fetched once from KEGG's own compound-name index, cached like everything
  else KEGG), not the bare C-number. Node SIZE is the node-importance measure
  chosen at setup: with relative betweenness centrality (the default), a
  bigger circle sits on the shortest reaction path between more OTHER
  compound pairs in this pathway - a hub or bottleneck the network routes
  through, while a small circle is structurally peripheral; with out-degree,
  a bigger circle directly feeds more reactions (it is a substrate more
  often). Both the NODES and their compound-name labels can be dragged with
  the mouse - move a node to untangle a crowded cluster (its edges and its
  own label follow it), or drag just a label when only the label needs to
  move. A dragged label's position survives switching tabs; "Reset label
  positions" puts every label (not moved nodes) back, and a "Compound
  labels" control hides them all if a large graph gets cluttered.

  Some compounds in the graph have no line to any other node. This is
  genuine KEGG data, not a rendering gap: KEGG's own diagram for that
  pathway places the compound (often a cofactor, a pathway endpoint, or
  something only linked via an enzyme/gene box the diagram does not encode
  as a reaction arrow) without connecting it to a reaction in THAT specific
  map. Verified live (2026-09-29): the human "Drug metabolism - other
  enzymes" map (hsa00983) has 63 compound nodes but only 30 participate in
  any reaction - 33 are genuinely isolated in KEGG's own figure. An isolated
  node still
  counts toward the node total and, if matched, toward Impact; it simply has
  zero centrality, like any dead end in the real network.


What to upload
--------------
Compound list (Over-Representation Analysis, Chemical Class, Lipid Class,
and Pathway Topology Analysis when you choose "compound list" there)

  A CSV or TSV with one compound per row.

    Required   Name        the compound name (e.g. Citrate, L-Glutamic acid,
                           PC 34:1)
    Optional   InChIKey, PubChem ID, KEGG ID, HMDB ID, ChEBI ID,
               LIPID MAPS ID, SMILES, Formula

  There are no statistics columns: no fold change, no p-value. A list has no
  "significant" subset to select - the compounds you upload ARE the test set.
  (Column names are matched case- and space-insensitively.)

  An optional BACKGROUND list (same format) narrows each database's reference
  set to compounds that appear in both the database and your background - for
  example the panel your instrument can actually measure.

Concentration table (Quantitative Enrichment, and Pathway Topology Analysis
when you choose "concentration table" there)

  Samples in rows: one sample-ID column, one group-label column (e.g. Group =
  Control / Treated), and every other column a compound's concentration or
  intensity. The transposed layout (compounds in rows, samples as columns, the
  group labels as a "Group" row) is also accepted.

  Both the layout and the group column are auto-detected and REPORTED back to
  you, with a control to override them - a wrong silent guess is worse than
  asking.

  Omiscera then computes, for each compound, a Welch t-test between the two
  groups you choose (Benjamini-Hochberg-adjusted), and the log2 fold change.
  You do not upload log2FC or p-values.
    - Log2 transform: on by default (concentration data are right-skewed).
      Zeros are handled by adding half the smallest positive value; the offset
      is reported. Untick it if your table is already log-scaled.
    - Missing values are NOT imputed. A compound needs at least the chosen
      number of measured samples (default 3) in EACH group to be tested;
      others are dropped and counted.
    - More than two groups: you choose which two to contrast. A multi-group
      test (ANOVA) produces an unsigned statistic that could not feed a
      direction-aware enrichment, so it is never used silently.
  The statistics table is shown on the results screen and downloadable, so the
  numbers Omiscera derived can be audited.

**Supply identifiers if you have them.** A compound name alone is genuinely
ambiguous: synonyms ("Palmitic acid" / "hexadecanoic acid" / "C16:0"),
isomers, and conjugate-base spellings ("Citrate" vs "Citric acid") all name
the same or nearly-the-same thing differently. Omiscera resolves names through
PubChem and KEGG, which handles most of this well, but an InChIKey or PubChem
CID you already have is always more reliable than a lookup.


The universe (what your list is tested AGAINST)
-----------------------------------------------
This is the most consequential design choice, and it differs by tile.

List tiles (ORA, Chemical Class, Lipid Class): the universe is the DATABASE's
own annotated set - every compound KEGG (or Reactome, RefMet, ...) annotates.
A plain list has no "measured" set, and using your own list as the universe
would put every one of your compounds in every term (a degenerate test). The
result tables say "universe: all N compounds annotated in <database>". If you
supply a background list, the universe becomes the compounds in both.

For KEGG the database's set is the organism's own maps (see KEGG below), and
Reactome/PathBank are likewise filtered to your organism.

Quantitative Enrichment: the universe is the compounds you actually measured.
Using a database's whole annotated set here would compare your terms against
compounds you never looked at.

Term size: only terms with at least the minimum (default 3) compounds in the
universe are tested, and, for the pathway tiles, at most the maximum (default
500 - a term holding a large share of the reference set says little). The two
class tiles have NO maximum by default: class taxonomies are big at exactly
the level people want (LIPID MAPS "Glycerophosphocholines" has thousands of
members), and a 500 cap silently emptied the very result the tile exists for.
Fisher's test is exact at any term size. Set a maximum on those tiles if you
want to drop the most generic classes (e.g. RefMet "Super class").


How a compound reaches each database
------------------------------------
This is the part worth understanding, because it explains most of what you
will see in the matching summary.

The databases are keyed on different identifiers, and two of them cannot be
reached from a name at all:

    KEGG pathways      KEGG C-number
    SMPDB              KEGG or ChEBI
    PathBank           KEGG or ChEBI
    Reactome           ChEBI only
    Chemical Class     RefMet (name/PubChem/InChIKey/KEGG/HMDB)
    Lipid Class        LIPID MAPS ID

So Omiscera builds a local cross-reference from bulk files (LIPID MAPS, RefMet
and ChEBI) and joins everything through one internal identifier. Your
compounds are matched into it in this order, most specific first:

    full InChIKey -> PubChem CID -> ChEBI ID -> KEGG ID -> HMDB ID
    -> LIPID MAPS ID -> name

and the matching summary reports how many compounds matched at each tier.
(Commas inside a compound name - "2,3-bisphosphoglycerate" - are handled: the
member lists in the result tables are joined with semicolons.)

Two tiers deserve attention:

  InChIKey skeleton (reported separately)
      An InChIKey has three blocks. The first encodes connectivity, the second
      stereochemistry and isotopes, the third protonation. When no exact match
      exists, Omiscera falls back to matching the first block alone. That
      bridges real and harmless differences - a measured lipid reported
      without sn-position, or an acid recorded as its anion - but it can also
      match a genuine stereoisomer of a different compound. These are counted
      and reported on their own so you can see how much of a result depends on
      them. (Example seen in testing: "Malate" matched the D- rather than the
      L- form this way.)

  Lipid shorthand (class-level only)
      A sum-composition name like "PC 34:1" names a bulk species covering many
      distinct structures. It cannot honestly be placed in a specific
      reaction, so these compounds are used for Chemical Class and Lipid Class
      enrichment and deliberately EXCLUDED from pathway membership. They are
      flagged class_level_only in the matched table. If you do know the exact
      structure, supply an InChIKey or LIPID MAPS ID and that row regains full
      pathway capability. In the Lipid Class analysis a shorthand compound is
      counted as one extra member of its own class (no single LIPID MAPS
      structure carries that name, so it could not otherwise appear at all).

This means a purely shorthand lipidomics upload will produce empty pathway
tables and full class tables. That is the correct result, not a failure - the
data supports the class question and not the pathway one.


Databases and licences
----------------------
  KEGG            Organism-specific: only the pathway maps YOUR organism has
                  are tested, reported under the organism's own ids and names
                  (hsa00020 "Citrate cycle (TCA cycle)", not the generic
                  map00020). By default KEGG's organism is derived from the
                  organism you chose (the same resolution Transcriptomics and
                  Proteomics use). A KEGG organism box on the setup screen -
                  a search over KEGG's own organism list, like the organism
                  box in Metabolic Connection - lets you pick one directly
                  instead, by name or code (eco, bsu, mmu...); an explicit
                  pick always wins and is safe to set apart from your main
                  organism, because Metabolomics has no gene-ID mapping that
                  could be resolved for one organism and tested against
                  another - the organism only selects pathway maps. When the
                  derivation is ambiguous (E. coli matches 64 KEGG genomes,
                  B. subtilis 14) KEGG is skipped and the box offers the
                  candidates. Numbers to expect, from KEGG itself: of 587
                  reference maps, human has 372, mouse 368, yeast 142,
                  Arabidopsis 162, E. coli 138 - a yeast run is no longer
                  tested against hundreds of pathways yeast does not have, and
                  the universe (N) shrinks accordingly, so p-values differ from
                  a generic-map run.
                  What the organism does and does NOT change. KEGG serves no
                  organism-specific compound lists (link/compound/hsa is
                  empty), and a map's compounds do not vary by organism
                  (verified: hsa00020 and sce00020 hold exactly the same 20
                  compounds as reference map00020). So the organism selects
                  WHICH MAPS are tested; the compounds inside a map are KEGG's
                  reference biochemistry. A pathway your organism has can
                  therefore contain a compound its enzymes never touch - a
                  small over-inclusion, stated in every KEGG result's caveat.
                  If KEGG has no genome for your organism, or the organism is
                  ambiguous and you have not confirmed one, KEGG returns an
                  EMPTY result with a note saying why. It never falls back to
                  the generic maps.
                  Pathway Topology Analysis only: unlike compound membership
                  above, a map's REACTIONS (which compound reacts with which)
                  ARE organism-specific - verified live: the human, yeast and
                  E. coli TCA cycle maps share the same 23 compounds but have
                  22, 21 and 21 reactions respectively. So the reaction graph
                  behind the Impact score is fetched and cached separately
                  per organism, even though the underlying compound list is
                  not.
                  KEGG's own map01000-map01099 range ("Overview of biosynthetic
                  pathways", "Biosynthesis of alkaloids derived from...", and
                  similar) is excluded: these are chemical-structure-class
                  meta-charts that aggregate compounds across huge swaths of
                  metabolism rather than naming a real pathway, and left in
                  they can outrank a genuinely specific hit (a test run saw
                  one claim EVERY one of its annotated compounds was
                  significant, because its real content was "most of central
                  metabolism"). Real, specific pathways that happen to
                  reference a metabolite only tangentially (a signaling
                  diagram citing a handful of small molecules as regulators)
                  are NOT excluded by this - a narrow experiment can still
                  produce a top hit that is technically correct but not the
                  most biologically central one; read the compound list on a
                  surprising top hit before trusting the ranking alone.
  Reactome        ChEBI2Reactome_All_Levels.txt, filtered by your organism.
                  Reactome and this module's own cross-reference frequently
                  cite DIFFERENT ChEBI IDs for the same compound at a
                  different protonation state (Reactome's CHEBI:16947
                  "citrate(3-)" versus this module's CHEBI:30769 "citric
                  acid" for the same uploaded "Citrate"), so exact-ID matching
                  alone misses most real matches. Bridged via the InChIKey
                  skeleton (see "Malate" above) for Reactome's own small,
                  closed set of referenced ChEBI IDs, fetched once through
                  ChEBI's REST API and cached; bridged compounds are counted
                  and reported in this table's own caveat. CC BY 4.0.
  SMPDB           Human only. Selecting another organism yields an explicit
                  empty result rather than a silent zero. CC BY 4.0.
                  SMPDB and PathBank contain thousands of drug-action and
                  disease pathways assembled from the same core metabolites, so
                  a list of central-carbon metabolites can put e.g. a
                  "<drug> Action Pathway" at the top. That is a statistically
                  honest overlap, not a claim about the drug - read the term
                  name, and prefer KEGG/Reactome for the core-metabolism
                  question.
  PathBank        Multi-organism. Uses the primary (canonical) pathway set by
                  default; the full set fragments the same biology across
                  hundreds of thousands of near-duplicate pathways, which
                  would swamp an FDR family. Open Database License (ODbL).
  Chemical Class  RefMet, from Metabolomics Workbench - a three-level taxonomy
                  (super class / main class / sub class) derived from
                  ClassyFire. **Metabolomics Workbench's terms of use restrict
                  it to personal, non-commercial use.** Omiscera ships code,
                  not data: the cache/ directory is gitignored and each
                  installation downloads for itself, so that restriction binds
                  your use of the results, not Omiscera's distribution. It
                  applies ONLY to the Chemical Class analysis - every other
                  source is freely licensed.
  Lipid Class     LIPID MAPS LMSD, using its own CATEGORY / MAIN_CLASS /
                  SUB_CLASS hierarchy. CC BY 4.0.
  ChEBI           Used only as a cross-reference source. CC BY 4.0.

If a source cannot be downloaded, the run continues with the others and says
which one was missing and what it cost. (ChEBI's bulk host, ftp.ebi.ac.uk, was
unreachable during development; the cross-reference still carried ~34,000
ChEBI IDs supplied by LIPID MAPS and RefMet, so Reactome matching survived.)

Please cite the databases you actually used, alongside Omiscera.


Reading the results
-------------------
Every result table is an independent multiple-testing family - each database
is corrected on its own, with no correction ACROSS databases (reporting the
best hit over all of them inflates the true error rate above any single FDR
column). A term is only tested if it contains at least one of your compounds
and clears the term-size limits; terms with no overlap are excluded entirely
rather than entered as p=1, because keeping them inflates the FDR for every
genuine test beside them.

List tiles - Fisher's exact test (one-sided, over-representation)
    k    your compounds inside the term
    K    all compounds in the term (in the universe)
    n    your compounds that are in the universe at all
    N    the universe (see above)
    Fold Enrichment, p-value, FDR, Bonferroni
    Overlapping Compounds - your own names, semicolon-separated

Quantitative Enrichment - competitive Wilcoxon rank-sum, two-sided
    Compounds in set, Mean score (set), Mean score (rest),
    Effect (mean difference), p-value, FDR, Bonferroni
    Each term's compounds are compared with every OTHER measured compound, by
    the chosen score (log2 fold change, or Welch t). There is no significance
    cutoff and no up/down split; the sign of the Effect gives the direction.
    This is NOT GSEA, deliberately: GSEA's permutation null assumes a long,
    dense ranked list (a transcriptome, ~20,000 genes), while a metabolomics
    panel is tens to a couple of thousand compounds whose missing values are
    not missing at random. The Proteomics module made the same choice.


Class tiles: what the figures show
----------------------------------
Chemical Class and Lipid Class results have a dot plot, a bar plot, a
compound-concept network, and a CLASS COMPOSITION pie.

The pie is descriptive: what fraction of your matched compounds falls in each
class, with a colour picker per class. It has no p-value - the statistical
result is the bar plot and the table. It is built from your whole resolved
compound set (not from the result table, which is filtered by term size and
significance and so would not sum to your compounds), counts each compound
once, and shows compounds with no class at the chosen level as "Unclassified"
rather than hiding them. "Classes shown" folds the smallest classes into one
"Other" slice.

Every slice gets a percentage label with a leader line to its slice. You can
hide the percentages, drag any label to move it (the line follows, and the
position is kept when the chart redraws), reset the positions, and use Zoom
(50-300%) to enlarge or shrink the whole figure. The labels are part of the
figure, so they appear in the SVG/PNG export where you left them.

A pie with many classes (dozens of RefMet sub classes) sizes itself: the figure
grows taller so every label has its own line and the whole legend is visible
without scrolling. "Label slices of at least" (0.5-5%) hides the labels of the
thinnest slices if you want a cleaner picture - they stay in the legend and in
the hover. "Slice borders" chooses white, black or no line between slices (the
hairline can dominate a pie with many thin slices).

Where the classes come from. Chemical Class uses RefMet's own taxonomy
(Metabolomics Workbench: super class / main class / sub class, derived from
ClassyFire and curated by RefMet) - Omiscera does not compute or guess a class.
Each compound of yours is matched to a RefMet record (by structure or ID,
falling back to name) and takes that record's class at the level you chose. A
compound RefMet does not contain has no class and shows as "Unclassified".
Lipid Class uses LIPID MAPS' own hierarchy the same way (Category / Main class
/ Sub class), and a lipid given only as shorthand ("PC 34:1") is placed by its
class prefix. A pie is only valid because every compound has exactly ONE
class at a given level, which is also why no equivalent is offered for the
pathway databases, where a compound sits in many pathways.

There is no enrichment map on these tiles. Class terms are disjoint (each
compound belongs to one class per level), so no two terms ever share a compound
and an enrichment map - which links terms that share compounds - could never
draw a single edge. The dot plot IS kept: it is a valid view of an enrichment
result for classes as much as for pathways.


Table layout (Quantitative Enrichment). Both common layouts are accepted and
detected automatically: samples in rows with one column per compound, or
compounds in rows with one column per sample. The group (condition) is a column
in the first layout and a row in the second, and it can be called anything
("Group", "Sample type", "Cohort", ...): the layout is recognised from the shape
of the table, and the Table layout / Group column selectors let you override
the guess. Blank rows and columns, compounds with no name and compounds with no
numeric value are ignored and counted; the same compound listed twice (for
example two adducts) is kept as two measurements and averaged into one score.
A table containing negative values cannot be log2-transformed - it is already
log-transformed or scaled - so untick "Log2-transform" for it.


Things worth keeping in mind
----------------------------
- Coverage is uneven. KEGG carries few individual lipid species, so a
  lipidomics dataset will match poorly to KEGG pathways and well to LIPID MAPS
  classes. This is a property of the databases, not of your data.
- Annotation bias applies here as it does everywhere: well-studied metabolites
  belong to more pathways, so they are easier to enrich for. A top hit can
  reflect study effort rather than biology.
- The class analyses test structural categories, not mechanisms. "Your
  compounds are mostly glycerophosphocholines" is a real and useful finding,
  but it is a statement about chemistry, not about a regulatory pathway.
- The first run downloads roughly 30 MB of reference files (and SMPDB's own
  file is large); later runs reuse the cache, which refreshes every 30 days.
