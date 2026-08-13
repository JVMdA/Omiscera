CUSTOM GENE-SET DATABASE — WHAT TO UPLOAD
===========================================

The ORA (Over-Representation Analysis) config screen has a "Custom database"
option. It lets you test your gene list against your own collection of gene
sets (e.g. an in-house pathway list) instead of, or alongside, the built-in
GO / Reactome / GWAS Catalog databases.

GENE IDENTIFIERS: UniProt ACCESSION *OR* Ensembl Gene ID
------------------------------------------------------------
The gene column in your file can hold EITHER:
  (a) UniProt accessions   e.g. P04637, Q00987, P38936
  (b) Ensembl Gene IDs      e.g. ENSG00000141510, ENSG00000012048

Name the gene column accordingly — "UniProt_Accession" for (a), or
"Ensembl_Gene_ID" for (b) — and select the matching option ("UniProt
accession" or "Ensembl Gene ID") in the upload widget on the ORA config
screen. Do not mix the two identifier types within one file.

The bundled example files below use UniProt accessions. To use Ensembl Gene
IDs instead, replace the values in the gene column with your Ensembl IDs
and rename that column to "Ensembl_Gene_ID" — no separate example file is
needed, just follow the same TIDY.csv / GMT.gmt layout with that column
swapped in.

If you pick the wrong option (e.g. your file has Ensembl IDs but you
selected "UniProt accession"), none of your genes will match and the
result shows a "0 genes matched" warning instead of a crash — so it's safe
to try, but re-upload with the correct selection to get real results.

WHAT IS A "GENE SET"?
----------------------
A gene set is just a named group of genes that belong together for some
reason — usually a pathway, a process, or a custom category you define.
Example: "p53_pathway" containing the genes TP53, MDM2, CDKN1A, ATM, ...

Your file is simply a list of these gene sets, each with:
  1. An ID/name for the set   (e.g. "p53_pathway")
  2. An optional description  (e.g. "p53 / DNA-damage response")
  3. Its member genes         (e.g. TP53, MDM2, CDKN1A, ...)

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
    gene / accession / ensembl / ensembl_gene_id (for the gene column) —
    but the names shown above are the clearest and recommended.
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
  example_custom_geneset_db_SYMBOLS.gmt  - GMT format, gene symbols (not
                                            Ensembl — shown for reference)

All three describe the SAME 10 example gene sets (p53/DNA-damage response,
cell cycle, apoptosis, PI3K-AKT signaling, MAPK/ERK signaling, Wnt
signaling, autophagy, DNA repair, immune/inflammatory response, and
glucose/insulin metabolism), built from real, verified human UniProt
entries — so they can be uploaded as-is to try the feature, or used as a
template for your own database.

NOTE: this README covers the UniProt accession / Ensembl Gene ID options.
If your uploaded gene list used Gene Symbols instead, select "Gene symbol"
in the upload widget and see README_custom_geneset_db.txt.
