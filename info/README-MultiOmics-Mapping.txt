Multi-omics — Cross-omics Pathway Mapping
=========================================

What this method is for
-----------------------
Joint Pathway Enrichment (see README-MultiOmics.txt) tells you WHICH pathways your
omics layers point at together. This method goes inside ONE pathway and shows
which of its members each omics hit, and how they connect: enzyme (gene / protein)
-> reaction -> metabolite. It is the "paint several omics on one pathway" view
that PaintOmics 4 (Liu et al. 2022), Pathview (Luo & Brouwer 2013) and KEGG Mapper
(Kanehisa & Sato 2020) popularised, built here on KEGG and Reactome.

IT IS DESCRIPTIVE. It computes no p-value. Whether a pathway is significant is
Joint Pathway Enrichment's question; this tells you what is going on inside it.

How to use it
-------------
1. Open the Multi-omics box "Cross-omics Pathway Mapping". Choose the omics (two
   to four), upload one file each, pick the organism, resolve the identifiers.
   If you already resolved your files in this session the box takes you straight
   on - nothing is uploaded twice. (The Joint Pathway Enrichment box always
   starts from a clean setup.)
2. Choose KEGG or Reactome and, for a differential-expression table, the padj and
   |log2 fold change| cutoffs. "List pathways hit by my omics" gives every pathway
   at least one omics hit, ranked by how many omics hit it and then by the total
   number of hits. "hit / annotated" is how many of the pathway's members that
   omics hit out of those it can place on the pathway. If you ran Joint Pathway
   Enrichment first, its combined FDR is shown beside each pathway. Search by name
   or ID, press Map.
3. You can also reach a pathway from a Joint Pathway Enrichment result: under its
   combined table, "Map this pathway".

What you get
------------
  Pathway map     The pathway drawn at the database's own positions and, for a
                  big map, at the database's own SIZE (see "Big maps"). A square
                  is a gene / protein (a KEGG box can hold several genes, a
                  Reactome box a whole complex or set), a circle is a metabolite.
  Members         Every member one of your omics measured or hit, one row per
                  omics: status, direction, log2FC, padj and HOW it was matched.
                  KEGG sometimes draws the same gene in several boxes; "Box ID"
                  tells them apart. A gene drawn in five boxes is ONE member hit.
  Reaction links  Reactions with evidence from more than one angle (below). Not
                  shown for a signalling map, which has no reactions.
  Interaction links  Only on a map that has relations (below).
  Overview        Counts, the descriptive caveat, and where a diagram came from.

READING THE MAP. Under each box sit up to four coloured chips, always in the same
order - Genomics, Transcriptomics, Proteomics, Metabolomics - so a box hit by
several omics never hides one:
    filled chip     the omics HIT this member
    hollow chip     the omics MEASURED it but it did not pass the cutoffs
    no chip         the omics did not have it (or it could not be placed here)
A triangle pointing up / down is the direction in a differential-expression table
(a diamond: the box's members disagree); a plain circle is a list hit (a plain
list has no direction, and no "measured but not significant" state). The colours
are the same as in the joint evidence plot. Hover a box or a chip for the members,
log2FC, padj and the identifier bridge. The spacing between chips is a fixed
number of screen pixels at every zoom level.

Three states, never two
-----------------------
"measured, not significant" is kept apart from "absent": a gene your RNA-seq saw
and that did not change is information, and a gene that was not in the table is
not. Only an omics that has a reference set - a differential-expression table -
can say "measured"; a plain list can only say hit or absent.

A box is hit if ANY of its members is. A KEGG box holding hsa:4967 and hsa:55753
is hit by Genomics when either gene is on your list.

How members are matched
-----------------------
Each match records which identifier bridge it went through (the "Matched by"
column and the hover).
  KEGG genes        the KEGG gene ID itself (Genomics, after the KEGG-ID step);
                    UniProt -> KEGG gene; Entrez Gene ID -> KEGG gene; and, as a
                    last resort, gene symbol -> KEGG gene. These are KEGG's own
                    conversion tables. An identifier that maps to MORE than one
                    KEGG gene is dropped, never guessed.
  KEGG compounds    the KEGG compound C-number from the Metabolomics matching.
  Reactome proteins UniProt directly (Genomics, Proteomics). A Transcriptomics
                    table has no UniProt column, so it goes Entrez -> KEGG gene
                    -> UniProt through KEGG's table when a KEGG organism can be
                    resolved for your organism, and by gene symbol only when it
                    cannot (and that is stated on the result).
  Reactome chemicals ChEBI ID first. Reactome often cites a different protonation
                    state of the same compound than a cross-reference does
                    (citrate(3-), CHEBI:16947, against citric acid, CHEBI:30769),
                    so a miss falls back to the InChIKey SKELETON - the same
                    bridge Metabolomics' own Reactome enrichment uses - and the
                    match is labelled as such.
A Reactome protein or complex is matched on its UniProt leaves only and a Reactome
chemical on its ChEBI leaves only, so a cofactor inside a complex (an iron-sulfur
cluster in ACO2) never makes the complex look like a metabolite hit.

Cross-omics reaction links
--------------------------
A reaction is listed when
  - its enzyme was hit by a gene / protein omics AND a substrate or product was
    hit by Metabolomics ("enzyme + metabolite"), or
  - its enzyme was hit in two or more omics ("enzyme in several omics").
Nothing else: an enzyme that was only measured, a metabolite whose enzyme nobody
hit, and a reaction with no evidence are not links. Direction is reported for the
expression omics only; the column says "concordant" or "discordant" only when
Transcriptomics AND Proteomics both hit the enzyme with a clear direction.
"Substrate" and "product" follow the database's nominal direction; "Reversible"
says when the reaction can run both ways. Descriptive - not a test.

Three kinds of map, and what they look like
-------------------------------------------
KEGG metabolic maps (Glycolysis, the TCA cycle, ...) are made of reactions: boxes
for enzymes, circles for compounds, lines for the reactions. Everything above
describes these.

KEGG's global overview maps ("Metabolic pathways" hsa01100 and the other 011xx /
012xx maps: 1,300+ genes, 3,500+ compounds) draw a gene not as a box but as a LINE
along its reaction arrow. Here that line IS the gene, and - as in KEGG's own
picture of the map - it is coloured by its pathway CLASS (Carbohydrate metabolism,
Lipid metabolism, Amino acid metabolism, Nucleotide metabolism, ...). KEGG stores that
colour in every line; its NAME comes from KEGG's own classification (br08901): a colour
is named by the metabolism class most of its genes' pathways belong to (the overview
maps themselves, and the disease pathways that share metabolic enzymes, do not vote).
What your omics did shows in the line's strength, and in the chips:
    faint line          no omics had this gene
    firmer line         measured, but not significant (a hollow chip)
    full colour, thick  HIT, with the same four chips as a box - Genomics,
                        Transcriptomics, Proteomics, Metabolomics, in fixed slots -
                        at the middle of the line
A pathway's label on the map is printed in its class's colour. Above the chart a legend
lists the classes (with how many of their lines were hit) and a colour picker for each:
choose another colour and exactly that class is redrawn; "KEGG colours" restores the
defaults. KEGG's own faint lines for enzymes this organism lacks are drawn as
decoration. So many labels would be unreadable: on a map with more than 300 boxes
only the hits are labelled; the "Labels" menu switches between hits only, hits and
measured, every box, and none. Hover a line (at its middle) for its members.

KEGG's signalling and disease maps (Apoptosis, NF-kappa B, the infection pathways,
...) have no reactions. Their lines are KEGG's RELATIONS between members:
    arrow, solid        activation, expression (and modifications that act on a
                        target: phosphorylation, ubiquitination, ...)
    bar                 inhibition, repression
    arrow, dashed       indirect effect
    plain line          binding / association, and compound-mediated relations
    dash-dot            dissociation
    dotted              missing interaction
    "+p" "-p" "+u"      phosphorylation, dephosphorylation, ubiquitination ...
A DARK relation is one whose BOTH ends were hit by some omics; a light one has at
least one end not hit. A group (a drawn complex) is an outline around its components
and is hit if ANY component is. The "Interaction links" pane lists the dark
relations: From, To, the relation, and which omics hit each end. It is KEGG's
statement of the relation plus your data's statement of who was hit - descriptive,
not a test.

Reactome pathways are drawn from Reactome's own diagram. Reactome's overview
pathways (Disease, Metabolism, Immune System, Infectious disease, Diseases of DNA
repair, ...: about one in five) have a diagram that holds ONLY links to their
sub-pathways, with no gene or compound on it. Mapping one gives no chart: it says so
and lists the sub-pathways with how many hits your omics have in each, ranked, with a
Map button on each, so one or two clicks reach a real diagram.

Big maps
--------
A map is drawn at the database's own scale - one screen pixel per KEGG / Reactome unit,
"100%" - and not squeezed into the width of the page. KEGG's global map is about
4,900 x 3,200 units: squeezed into a 1,000-pixel page its circles, lines and labels
land on top of each other (measured: three quarters of the hit labels sat on another
label), at 100% about a tenth do. A map wider than the page is therefore drawn at 100%
inside a frame you scroll with the scroll bars or the wheel; a smaller one is fitted to
the page, as before. The "Size" menu changes it: fit to page, 50%, 75%, 100%, 150%,
200%. The chart's own toolbar (zoom, pan, autoscale, reset axes, download) sits at the
chart's right edge, which is far off-screen on a 4,900-pixel chart. So the same tools are
repeated in a row above it whenever the map is drawn at a fixed size: Zoom tool (drag a
box to zoom to it), Pan tool, Zoom in, Zoom out, Autoscale (fit all the data), Reset axes
(back to the view the map opened with) and "Download map" (SVG, PNG or TIFF). When the
map is fitted to the page the chart's own toolbar is in view and the row is hidden.

Where the pathway structure comes from
--------------------------------------
KEGG       the pathway's KGML file: its genes, compounds, reactions, relations
           and groups with KEGG's own layout (rest.kegg.jp/get/<pathway>/kgml),
           including the line geometry of the global maps. Orthologs are drawn but
           carry no members; a glycan is never a compound hit.
Reactome   the pathway's diagram layout plus its entity graph (a complex or set is
           opened down to its UniProt and ChEBI leaves). A pathway with no
           diagram of its own - Reactome draws many sub-pathways inside a parent -
           is drawn from its NEAREST ancestor that has one, restricted to the
           pathway's own reactions, and the Overview says which ancestor. An
           overview pathway (above) is detected from its diagram and its
           sub-pathways are looked up one by one.
Both are cached on disk after the first fetch. A transport error is reported and
never cached; "this pathway has no map" is cached.

Why not MetaCyc
---------------
MetaCyc was on the original wish list. Since 1 January 2024 it requires a paid
BioCyc subscription for the website, the downloads and the web-services API
(biocyc.org/metacyc-paid.shtml). A free, public Omiscera cannot depend on it, so
it is not offered. KEGG and Reactome cover genes AND compounds, which is what a
map of enzymes and metabolites needs.

Limits to keep in mind
----------------------
  - A hit is only as good as the cutoffs. The map shows what passed them; change
    them and list the pathways again.
  - Matching is by identifier, not by name. A compound without a KEGG / ChEBI /
    InChIKey match, and a gene KEGG's tables cannot place, do not appear (the
    resolve screen counts what was left out).
  - A KEGG map only holds the genes KEGG assigns to that reaction in your
    organism; a Reactome diagram of an organism other than human is Reactome's
    own inference from human.
  - Large top-level pathways make large maps; the biggest (hsa01100, about 5,000
    boxes) is a few thousand pixels across, draws in under a second once its data is
    here (the first request for it also downloads KEGG's file), and is meant to be
    scrolled around. Labels can still touch in the densest corners; the Labels menu
    and a larger Size help. There is no dragging of boxes.
  - A relation in a KEGG signalling map is drawn between the two boxes it names; where
    KEGG draws the same gene in several boxes, each box has its own relations.
  - A reaction link is evidence to look at, not a finding: an enzyme and its
    product both being hit does not by itself say the flux changed.
  - The picker is a coverage list, not a test; its ranking says where most of
    your data lands, not what is significant.

Downloads
---------
The member table and the reaction-link table of the pathway on screen download as
CSV with every row (the on-screen tables only scroll).

Citations
---------
  Liu T, Salguero P, Petek M, et al. (2022). PaintOmics 4: new tools for the
      integrative analysis of multi-omics datasets supported by multiple pathway
      databases. Nucleic Acids Research 50(W1):W551-W559.
  Luo W, Brouwer C (2013). Pathview: an R/Bioconductor package for pathway-based
      data integration and visualization. Bioinformatics 29(14):1830-1831.
  Kanehisa M, Sato Y (2020). KEGG Mapper for inferring cellular functions from
      protein sequences. Protein Science 29(1):28-35.
  Fabregat A, Sidiropoulos K, Viteri G, et al. (2018). Reactome diagram viewer:
      data structures and strategies to boost performance. Bioinformatics
      34(7):1208-1214.
  Milacic M, Beavers D, Conley P, et al. (2024). The Reactome Pathway
      Knowledgebase 2024. Nucleic Acids Research 52(D1):D672-D678.
  Kanehisa M, Furumichi M, Sato Y, et al. (2023). KEGG for taxonomy-based
      analysis of pathways and genomes. Nucleic Acids Research 51(D1):D587-D592.
