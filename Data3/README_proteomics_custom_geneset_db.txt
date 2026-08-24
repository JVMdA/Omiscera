PROTEOMICS CUSTOM GENE-SET DATABASE — WHAT TO UPLOAD
=======================================================

The Proteomics Over-Representation Analysis config screen has a "Custom
database" option. It lets you test your significant proteins against your
own collection of gene sets (e.g. an in-house pathway list) instead of, or
alongside, the built-in GO / Reactome / MSigDB / InterPro / Complex Portal /
UniProt databases.

GENE IDENTIFIERS: UniProt ACCESSION ONLY
-------------------------------------------
Unlike the Genomics module's custom database, there is no gene-symbol or
Ensembl option here. The gene column in your file must hold UniProt
accessions, e.g. P04637, Q00987, P38936 - the same identifier space every
Proteomics analysis uses throughout (see the main Proteomics README for
why: MS search engines emit UniProt accessions natively, so there is
nothing to "resolve" the way a typed gene name needs to be).

If your file has gene symbols instead, it will still be accepted and
parsed - your custom database "loaded" with real term/gene counts - but
none of those symbols will match any of your resolved protein accessions,
so the run will show "0 genes matched" instead of a crash. Convert your
gene column to UniProt accessions first (e.g. via UniProt's own ID mapping
tool at uniprot.org/id-mapping) if you see that warning.

WHAT IS A "GENE SET"?
------------------------
A gene set is just a named group of proteins that belong together for some
reason - usually a complex, a pathway, or a custom category you define.
Example: "26S_proteasome" containing the 33 real subunits of the 26S
proteasome.

Your file is simply a list of these gene sets, each with:
  1. An ID/name for the set     (e.g. "26S_proteasome")
  2. An optional description    (e.g. "26S proteasome complex")
  3. Its member accessions      (e.g. O00231, O00232, O00487, ...)

TWO ACCEPTED FILE FORMATS
----------------------------

FORMAT 1 - Tidy table (CSV or TSV), one gene per row
See: example_proteomics_custom_geneset_db_TIDY.csv

    Gene_Set_ID,Gene_Set_Name,UniProt_Accession
    26S_proteasome,26S proteasome complex,O00231
    26S_proteasome,26S proteasome complex,O00232
    20S_proteasome_core,20S proteasome complex,O14818
    ...

  - Column names can also be: term_id / pathway / gene_set (for the ID),
    term_name / description / pathway_name (for the name), and
    gene / accession / uniprot_id (for the gene column) - but the names
    shown above are the clearest and recommended.
  - One row per (gene set, gene) pair. A gene set with 33 members needs 33
    rows.

FORMAT 2 - GMT (Gene Matrix Transposed), one gene set per line
See: example_proteomics_custom_geneset_db.gmt

    Gene_Set_ID <TAB> Gene_Set_Name <TAB> Gene1 <TAB> Gene2 <TAB> Gene3 ...

    26S_proteasome	26S proteasome complex	O00231	O00232	O00487	...

  - This is the standard format used by GSEA / MSigDB / Enrichr.
  - Columns are separated by TABS, not commas. Each line is one gene set,
    with as many genes as needed after the first two columns.

Both formats describe exactly the same information - pick whichever is
easier to build (a spreadsheet naturally exports the tidy CSV; GMT is more
compact and is the standard exchange format if you already have one).

BACKGROUND / UNIVERSE
------------------------
Your custom database's own genes are still narrowed to the proteins
actually QUANTIFIED in your experiment before testing (same convention as
every other database here) - the statistical background is never simply
"everything in your file," it is the intersection of your file's genes
with your resolved protein table.

EXAMPLE FILES IN THIS FOLDER
--------------------------------
  example_proteomics_custom_geneset_db_TIDY.csv  - tidy table
  example_proteomics_custom_geneset_db.gmt       - GMT format

Both describe the SAME 6 example gene sets - the 26S and 20S proteasome
complexes, the 40S and 60S cytosolic ribosomal subunits, and two
spliceosomal complexes (Pre-B and B) - built from real Complex Portal
entries, the same six complexes seeded into the bundled
example_proteomics_de_full.csv DE example (each complex was given its own,
different rate of being significant, matching how real regulation rarely
hits every complex equally). Uploading this custom database against that
example is a useful way to sanity-check that your own custom database is
being parsed and matched correctly before trusting it on real data - note
that a custom database's statistical background is the union of ONLY the
genes inside your own file, so it tests these six complexes AGAINST EACH
OTHER, not against the whole proteome. Only the complex(es) with the
highest designed hit-rate (the ribosomal subunits) are expected to come
back significant here - that is the correct, informative result for this
kind of comparison, not a sign something is missing. To see every one of
these six complexes individually confirmed as enriched against the WHOLE
proteome, look at the built-in Complex Portal database's own results on
the same example instead.
