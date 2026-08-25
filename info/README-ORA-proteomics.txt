PROTEOMICS — OVER-REPRESENTATION ANALYSIS (ORA)
==================================================

WHAT THIS MODULE IS FOR
--------------------------
You give Omiscera a protein-level differential-abundance table — every
protein you quantified, with its fold-change and p-value. ORA first applies
your padj/log2FC thresholds to decide which proteins count as significantly
changed, then asks:

    "Is there a common biological theme among my significant proteins — do
    they belong to a specific process, pathway, complex, domain family or
    disease more than I'd expect if I'd picked proteins at random from
    everything I actually quantified?"

This is repeated for every term in each selected database, separately for
all significant proteins, up-regulated ones, and down-regulated ones,
producing a ranked table of over-represented biological themes per
direction.

WHAT YOU UPLOAD
------------------
A protein-level differential-abundance table: a UniProt accession column, a
log2 fold-change column, and a p-value column. An adjusted p-value
(padj/FDR) column is optional — it is computed for you via
Benjamini-Hochberg if absent.

This is what you already have coming out of Perseus, limma, MSstats, or any
equivalent differential analysis. You do NOT need to look up UniProt
accessions: your MS search engine searched a protein FASTA (almost always
UniProt-derived) and already wrote accessions into its output —
MaxQuant's "Protein IDs"/"Majority protein IDs", DIA-NN's "Protein.Group",
FragPipe's "Protein ID", and so on.

Three things are handled for you, and none of them have a transcriptomics
equivalent:

  - PROTEIN GROUPS. MS often cannot distinguish isoforms or close homologs,
    so tools report a group like "P04637;P04637-2;Q9H3D4". The leading
    accession is taken as the group's representative, which is the standard
    convention.
  - ISOFORM SUFFIXES. "P04637-2" is reduced to "P04637", because annotation
    databases key on the canonical entry.
  - CONTAMINANTS AND DECOYS. Rows whose identifier carries MaxQuant's
    CON__/REV__ prefix, AND rows flagged in a "Potential contaminant",
    "Reverse" or "Only identified by site" column, are removed BEFORE
    anything else, and the counts are reported on screen. Left in, these
    silently corrupt every downstream result — a contaminant is never a
    real biological finding. This matches what Perseus does.

Deprecated or merged UniProt accessions are also resolved to their current
primary accession (UniProt periodically retires entries), and two rows that
resolve to the same protein are collapsed, keeping the stronger evidence.

Identifier space is UniProt ACCESSION only — there is no gene-symbol option,
unlike Genomics. A symbol would collapse the isoform distinction MS just
made, and a protein group can span several genes.

THE UNIVERSE (BACKGROUND)
------------------------------
The universe is every protein actually QUANTIFIED in this experiment, not a
generic whole-proteome background — then narrowed further, per database, to
those also annotated in that particular database. A protein that mass
spectrometry never detected could not have been called significant
regardless of its biology, so it does not belong in the background either.

This matters more in proteomics than in transcriptomics: MS detects only a
subset of the proteome (typically 3,000-8,000 proteins), so using a
whole-proteome background would treat thousands of never-measured proteins
as if they had been measured and found unchanged. That inflates every
p-value dramatically — measured on this app's own bundled example, using the
whole GO:CC-annotated proteome (~150,000 entries including isoforms and
unreviewed entries) instead of the ~3,300 quantified-and-annotated proteins
made the top hit's p-value roughly 10^33 times smaller, and even changed
which term ranked first.

An optional uploaded background list narrows the universe further, e.g. to a
targeted panel. It applies globally to every selected database.

THE STATISTICAL TEST
-----------------------
A one-sided Fisher's exact test (hypergeometric), the same core test as
Genomics and Transcriptomics ORA.

No gene-length bias correction is offered here. That correction exists in
Transcriptomics because longer transcripts accumulate more reads and are
more likely to be called DE. The dominant detection bias in MS is protein
ABUNDANCE, not length, and there is no established proteomics equivalent of
goseq's correction — so rather than apply a transcriptomics correction that
is only approximately right, none is applied, and this is stated rather than
left for you to assume.

Terms that contain none of your significant proteins are never tested at
all. They cannot produce a finding (their p-value is 1.0 by construction),
and including them would only inflate the multiple-testing family and
worsen every other term's FDR. This follows the clusterProfiler convention.

FOLD ENRICHMENT
------------------
Alongside the p-value, each term gets a fold enrichment: how many times more
often proteins from this term appear among your significant proteins than
chance (given the universe) would predict. Fold enrichment measures how
STRONG the enrichment is; the p-value measures how CONFIDENT you can be it
isn't a coincidence — complementary, not interchangeable.

THE DATABASES
----------------
Eight sources are offered. Six are natively keyed on UniProt accessions —
they annotate PROTEINS, not genes — so your accessions are used directly,
with no conversion:

  - GENE ONTOLOGY (GO) — what a protein does, where it acts, and what
    process it takes part in. Tested as three separate families (Biological
    Process, Molecular Function, Cellular Component). Source: EBI QuickGO.
  - REACTOME — curated metabolic and signaling pathways.
  - INTERPRO — protein domain and family signatures. Genuinely
    protein-level: a domain is a property of a protein sequence. Pfam,
    SMART, PROSITE and others are member databases INSIDE InterPro, so they
    are already covered by it.
  - COMPLEX PORTAL (EMBL-EBI, CC0) — known macromolecular complexes.
    Proteins physically act in complexes, which makes this one of the most
    proteomics-native questions you can ask. Covers ~20 well-studied
    species.
  - UNIPROT (DISEASE) — UniProt's own curated disease annotations.
  - UNIPROT (KEYWORDS) — UniProt's curated functional/technical keywords.

Two are gene-set resources keyed on gene SYMBOLS, so they are tested through
the gene symbol fetched for each of your proteins during resolution:

  - MSIGDB — 8 curated collections you pick from (human and mostly mouse).
    See the Transcriptomics ORA README for the full per-collection
    availability and licensing table.
  - DISEASES (Jensen Lab, CC BY 4.0, HUMAN ONLY) — curated gene-disease
    associations. See the next section for why this exists alongside
    UniProt's own disease annotations.

Plus your own uploaded sets:

  - CUSTOM DATABASE — GMT or a tidy term/gene table, in UniProt-accession
    space. An example (both formats, with a README) is downloadable from the
    config screen. A file in gene symbols is rejected up front with an
    explanation rather than silently matching nothing.

NOT INCLUDED: PANTHER. Its website and FTP server sit behind bot protection
that blocks automated download entirely (every request returns 403 or a
challenge page, including for static files). This is an access barrier, not
a licensing one. PANTHER's protein FAMILY signatures are in any case already
covered through InterPro, which includes PANTHER as a member database; only
PANTHER's own Protein Class hierarchy is missing.

WHY TWO DISEASE DATABASES
-----------------------------
They are not redundant, and the difference is worth understanding before
reading either one's results.

UniProt's disease annotations are overwhelmingly MONOGENIC — one gene, one
named rare disease. Measured on the human set: 97.5% of its 6,862 disease
terms contain exactly ONE gene ("Diamond-Blackfan anemia 22",
"Diamond-Blackfan anemia 19", and so on as separate entries). Only 77 terms
(1.1%) contain 3 or more. A one-gene term cannot really be "enriched" — it
is just restating that one protein was significant — so nearly the whole
database is correctly filtered out by the minimum-term-size setting. That is
not a defect; it is what the resource is.

DISEASES aggregates comparable knowledge into disease CATEGORIES instead:
396 of its 1,115 curated terms (35.5%) contain 3 or more genes, and its
largest terms are things like "Neurodegenerative disease" (442),
"Cancer" (326) and "Epilepsy" (324). That is what makes it usable for
enrichment where UniProt's disease entries are not.

Both are kept, because UniProt's narrow entries are still exactly right if
you are studying one specific rare disease. Only the manually-curated
channel of DISEASES is used; its text-mined and GWAS-scale experimental
channels were deliberately left out as too low-confidence and too broad
respectively.

CORRECTING FOR TESTING MANY TERMS AT ONCE
---------------------------------------------
Every term tested is an opportunity for a false positive, so each result
table carries both a Benjamini-Hochberg FDR and the stricter Bonferroni
value alongside the raw p-value.

Each database x direction is its own INDEPENDENT family — GO Biological
Process (up) is corrected separately from Reactome (all), and each selected
MSigDB collection is its own family too. There is no correction ACROSS
families. This means scanning every table and reporting whichever hit looks
best inflates your true error rate above what any single FDR column states.
The results screen says how many independent families a run produced.

RESULTS AND VISUALIZATION
-----------------------------
  - Summary cards — universe size and the number of significant proteins,
    all/up/down.
  - Result table — one row per tested term: Term ID, Term Name, k, K, n, N,
    fold enrichment, p-value, FDR, Bonferroni, and the overlapping
    proteins. Downloadable as CSV with full numeric precision.
  - "Overlapping Proteins" lists each member as ACCESSION (SYMBOL), e.g.
    "P05386 (RPLP1)" — the accession is the unambiguous identifier, the
    symbol is what makes it readable. Every database uses this same format,
    including the two symbol-keyed ones.
  - Figures — dot plot, bar plot, protein volcano, enrichment map, and a
    protein-concept network coloured by log2FC. All exportable.

  NOTE: the figures show only terms passing your significance threshold,
  while the result TABLE lists every term that was tested. This is
  deliberate — the table is the full statistical record, the figures are the
  findings — but it does mean a database can legitimately show, say, 26
  tested terms and only 4 plotted points. Raise the "Significance threshold
  (FDR <)" box to see more.

THINGS TO KEEP IN MIND
--------------------------
  - ORA only ever looks at proteins that passed your cutoffs. A coordinated
    but individually sub-threshold shift across a whole pathway is invisible
    to it.
  - Your universe is what MS actually detected, which is biased toward
    abundant proteins. Low-abundance regulatory proteins (kinases,
    transcription factors, receptors) are systematically under-represented
    in any shotgun proteomics experiment, and therefore in these results.
  - The two symbol-keyed databases (MSigDB, DISEASES) cannot distinguish
    protein isoforms — one gene is one entity to them. If your data
    separated two isoforms of the same gene, that distinction is lost for
    those two databases only.
  - DISEASES is human-only and MSigDB is human/mouse-only; both are greyed
    out for other organisms rather than returning an unexplained empty
    table. Complex Portal covers ~20 species. GO, Reactome, InterPro and
    UniProt cover far more.
  - A database returning no terms is a real answer, not necessarily a
    problem — check how many terms were TESTED (shown in the selector)
    before assuming something went wrong.
