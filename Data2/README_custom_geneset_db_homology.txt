CUSTOM GENE-SET DATABASE — WHAT TO UPLOAD
===========================================

The ORA (Over-Representation Analysis) config screen has a "Custom database"
option. It lets you test your gene list against your own collection of gene
sets (e.g. an in-house pathway list) instead of, or alongside, the built-in
GO / Reactome / GWAS Catalog databases.

GENE IDENTIFIERS: UniProt ACCESSION ONLY (Homology / BLAST input)
------------------------------------------------------------------
You uploaded your gene list as SEQUENCES (Homology / BLAST). BLAST resolves
each sequence to a UniProt accession, so a custom database used with this
input must identify its genes by UniProt accession — e.g. P04637, Q00987,
P38936.

Name the gene column "UniProt_Accession" and leave the id-type selector on
"UniProt accession" (it is the only valid option for BLAST input). Gene
symbols and Ensembl Gene IDs are NOT available here, because BLAST does not
produce them.

The bundled example files below already use UniProt accessions — you can
upload them as-is to try the feature, or use them as a template for your own
database.

WHAT IS A "GENE SET"?
----------------------
A gene set is just a named group of genes that belong together for some
reason — usually a pathway, a process, or a custom category you define.
Example: "p53_pathway" containing the genes TP53, MDM2, CDKN1A, ATM, ...

Your file is simply a list of these gene sets, each with:
  1. An ID/name for the set   (e.g. "p53_pathway")
  2. An optional description  (e.g. "p53 / DNA-damage response")
  3. Its member genes         (as UniProt accessions)

TWO ACCEPTED FILE FORMATS
--------------------------

FORMAT 1 — Tidy table (CSV or TSV), one gene per row
See: example_custom_geneset_db_TIDY.csv

    Gene_Set_ID,Gene_Set_Name,UniProt_Accession
    p53_pathway,p53 / DNA-damage response,P04637
    p53_pathway,p53 / DNA-damage response,Q00987
    cell_cycle,Cell cycle progression,P24385
    ...

  - Column names can also be: term_id / pathway / gene_set (for the ID),
    term_name / description / pathway_name (for the name), and
    gene / accession / uniprot_id (for the gene column) — but the names
    shown above are the clearest and recommended.
  - One row per (gene set, gene) pair. A gene set with 9 genes needs 9 rows.

FORMAT 2 — GMT (Gene Matrix Transposed), one gene set per line
See: example_custom_geneset_db_GMT.gmt

    Gene_Set_ID <TAB> Gene_Set_Name <TAB> Gene1 <TAB> Gene2 <TAB> Gene3 ...

    p53_pathway	p53 / DNA-damage response	P04637	Q00987	P38936	...
    cell_cycle	Cell cycle progression	P24385	P11802	Q00534	...

  - This is the standard format used by GSEA / MSigDB / Enrichr.
  - Columns are separated by TABS, not commas. Each line is one gene set,
    with as many genes as needed after the first two columns.

Both formats describe exactly the same information — pick whichever is
easier to build (a spreadsheet naturally exports the tidy CSV; GMT is more
compact and is the standard exchange format if you already have one).

BACKGROUND / UNIVERSE
-----------------------
No species or background list is needed for a custom database — the
statistical background (the "universe") is automatically the union of
every gene that appears anywhere in your file.

EXAMPLE FILES IN THIS FOLDER
-------------------------------
  example_custom_geneset_db_TIDY.csv     - tidy table, UniProt accessions
  example_custom_geneset_db_GMT.gmt      - GMT format, UniProt accessions

Both describe the SAME 10 example gene sets (p53/DNA-damage response, cell
cycle, apoptosis, PI3K-AKT signaling, MAPK/ERK signaling, Wnt signaling,
autophagy, DNA repair, immune/inflammatory response, and glucose/insulin
metabolism), built from real, verified human UniProt entries.
