# -*- coding: utf-8 -*-
"""
Genomics — Gene Function Annotation (GO + UniProt).

Unlike ORA / Burden Test / Network Propagation (all group-level: "what pattern
emerges from this gene list as a whole?"), this is a per-gene DESCRIPTIVE
lookup: for each resolved gene, return its Gene Ontology terms (Biological
Process / Molecular Function / Cellular Component) and its UniProt functional
description, protein name, and keywords. No p-value, no enrichment ratio, no
background universe — nothing is tested against chance, so there is no
statistics engine here, only fetch + classify + tabulate.

Data sources (native, not the bulk downloads ORA/Reactome/GWAS use — this is a
targeted PER-PROTEIN lookup, not a whole-organism universe build):
  - GO   : QuickGO annotation/search REST API (EBI), queried by UniProt
           accession (geneProductId), which returns every GO annotation for
           that protein with its evidence code — taxon-independent, so it also
           powers single-gene UniProt-ID mode with no organism selection.
  - UniProt: uniprotkb/search REST API, requesting the structured FUNCTION
           comment (already evidence-tag-free prose, unlike the flat-file
           format), keywords, and review status.

Differentiator: `classify_annotation_confidence` is a per-gene "how
well-characterized is this gene, really?" indicator (well-characterized /
moderate / sparse), built from GO evidence quality + UniProt curation. GO
annotation is famously skewed — the majority of annotations cover a minority
of genes (Wood et al. 2012, PLOS Comput Biol) — so any gene list is a mix of
well-studied and "dark" genes, and this makes that visible at a glance instead
of requiring the user to manually judge each gene's evidence codes.
"""

import re
import time

import pandas as pd

import Genomics_module as genomics_ora  # HTTP_SESSION / OMISCERA_HEADERS reuse

# ---------------------------------------------------------------------------
# Single-gene mode: identifier-type auto-detection
# ---------------------------------------------------------------------------

# Canonical UniProt accession pattern (6-char: [OPQ] prefix; 6- or 10-char:
# [A-NR-Z] prefix — O/P/Q are reserved for the first form). Matching this means
# resolution can be skipped entirely, which is also what makes single-gene
# mode work with no organism selected.
UNIPROT_ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
# Ensembl gene ID: ENS + optional species infix (mouse ENSMUSG, etc.) + G + 11 digits.
ENSEMBL_GENE_ID_RE = re.compile(r"^ENS[A-Z]*G[0-9]{11}(\.[0-9]+)?$")


def looks_like_uniprot_accession(identifier: str) -> bool:
    return bool(UNIPROT_ACCESSION_RE.match(str(identifier).strip().upper()))


def looks_like_ensembl_gene_id(identifier: str) -> bool:
    return bool(ENSEMBL_GENE_ID_RE.match(str(identifier).strip().upper()))


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

QUICKGO_ANNOTATION_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
QUICKGO_PAGE_LIMIT = 200          # verified live: QuickGO's max page size (limit>200 -> HTTP 400)
QUICKGO_MAX_PAGES_PER_BATCH = 100  # safety cap: 100 * 200 = 20,000 annotations/batch, generous
QUICKGO_BATCH_SIZE = 15           # accessions per geneProductId= comma-list

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FUNCTIONAL_FIELDS = (
    "accession,protein_name,organism_name,cc_function,keyword,reviewed,protein_existence"
)
UNIPROT_BATCH_SIZE = 25

DEFAULT_TIMEOUT = 30
DEFAULT_RETRIES = 3

# GO evidence codes (short "goEvidence" form, e.g. "IDA" — NOT the ECO:XXXXXXX
# id, which QuickGO also returns but is not human-scannable). Grouped per the
# GO Consortium's own evidence-code categories (geneontology.org/docs/guide-to-go-evidence-codes).
EXPERIMENTAL_EVIDENCE_CODES = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP",
}
CURATED_NONEXPERIMENTAL_EVIDENCE_CODES = {
    "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "TAS", "NAS", "IC",
}
# Not filtered on here (see run_gene_function_annotation's docstring) - this is
# the code value the frontend's live "Include IEA" toggle filters Evidence
# Code == "IEA" against, kept here as the single documented source of truth.
IEA_EVIDENCE_CODE = "IEA"

ANNOTATION_COLLAPSED_COLUMNS = [
    "Gene", "UniProt ID", "Protein Name", "Organism", "Reviewed",
    "GO Biological Process", "GO Molecular Function", "GO Cellular Component",
    "Keywords", "UniProt Function", "Annotation Status", "Experimental GO Terms", "Total GO Terms",
]
ANNOTATION_TIDY_COLUMNS = [
    "Gene", "UniProt ID", "GO ID", "Term", "Namespace", "Evidence Code", "Qualifier",
]

GO_ASPECT_TO_COLUMN = {
    "biological_process": "GO Biological Process",
    "molecular_function": "GO Molecular Function",
    "cellular_component": "GO Cellular Component",
}

STATUS_LABELS = {
    "well_characterized": "Well-characterized",
    "moderate": "Moderate",
    "sparse": "Sparse (dark)",
}

# In-process per-accession caches (mirrors other retriever/engine caches in
# this codebase) so re-running with a different IEA toggle, or annotating an
# overlapping gene list, never re-fetches an accession already seen this run.
_GO_CACHE: dict = {}
_UNIPROT_FUNCTIONAL_CACHE: dict = {}


# ---------------------------------------------------------------------------
# GO annotations (QuickGO, per-protein)
# ---------------------------------------------------------------------------

def _quickgo_annotation_get(params: dict, verbose: bool = False) -> "dict | None":
    """Single GET against QUICKGO_ANNOTATION_URL with retry. Returns parsed
    JSON dict or None on failure."""
    if genomics_ora.HTTP_SESSION is None:
        return None
    headers = dict(genomics_ora.OMISCERA_HEADERS)
    headers["Accept"] = "application/json"
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            resp = genomics_ora.HTTP_SESSION.get(
                QUICKGO_ANNOTATION_URL, params=params, headers=headers, timeout=DEFAULT_TIMEOUT,
            )
            if resp.status_code == 200:
                return resp.json()
            last_error = f"HTTP {resp.status_code}"
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < DEFAULT_RETRIES:
                time.sleep(1.5)
                continue
            break
        except Exception as exc:
            last_error = str(exc)
            if attempt < DEFAULT_RETRIES:
                time.sleep(1.0)
    if verbose:
        print(f"[GeneFunctionAnnotation] QuickGO request failed: {last_error}")
    return None


def fetch_go_annotations(
    accessions: list,
    batch_size: int = QUICKGO_BATCH_SIZE,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Fetch EVERY GO annotation (all evidence, including IEA) for each
    UniProt accession via QuickGO's per-protein annotation/search endpoint.

    Deliberately NOT the taxon-wide downloadSearch universe builder used by
    ORA (Genomics_module.build_go_universe): that fetches the whole
    organism's annotations, which is wasteful for a handful of genes and
    cannot work at all in single-gene UniProt-ID mode (no taxon).

    All evidence is fetched (never pre-filtered by IEA) so a display-only
    "include IEA" toggle can filter downstream with no re-query, and so
    classify_annotation_confidence always sees the full evidence picture.

    Returns {accession: [{"go_id", "name", "namespace", "evidence", "qualifier"}, ...]}.
    Negated annotations (qualifier starting with "NOT") are excluded — a
    "NOT|enables" annotation means the gene does NOT have that function, so
    including it as a positive annotation would be actively wrong.
    """
    result: dict = {acc: _GO_CACHE[acc] for acc in accessions if acc in _GO_CACHE}
    to_fetch = [acc for acc in dict.fromkeys(accessions) if acc not in _GO_CACHE]
    total = len(to_fetch)
    done = 0

    for i in range(0, len(to_fetch), batch_size):
        batch = to_fetch[i:i + batch_size]
        batch_set = set(batch)
        for acc in batch:
            result.setdefault(acc, [])

        page = 1
        total_pages = 1
        while page <= total_pages and page <= QUICKGO_MAX_PAGES_PER_BATCH:
            params = {
                "geneProductId": ",".join(batch),
                "limit": QUICKGO_PAGE_LIMIT,
                "page": page,
                "includeFields": "goName",
            }
            data = _quickgo_annotation_get(params, verbose=verbose)
            if data is None:
                break
            page_info = data.get("pageInfo") or {}
            total_pages = int(page_info.get("total") or 1)
            for row in data.get("results") or []:
                gp_id = str(row.get("geneProductId", ""))
                acc = gp_id.split(":", 1)[1] if ":" in gp_id else gp_id
                if acc not in batch_set:
                    continue
                qualifier = str(row.get("qualifier") or "")
                if qualifier.upper().startswith("NOT"):
                    continue
                result[acc].append({
                    "go_id": str(row.get("goId", "")),
                    "name": str(row.get("goName") or ""),
                    "namespace": str(row.get("goAspect", "")),
                    "evidence": str(row.get("goEvidence", "")),
                    "qualifier": qualifier,
                })
            page += 1

        for acc in batch:
            _GO_CACHE[acc] = result.get(acc, [])
        done += len(batch)
        if progress_callback:
            try:
                progress_callback(done, total, f"Fetching GO annotations ({done}/{total} genes)")
            except Exception:
                pass

    return {acc: result.get(acc, []) for acc in accessions}


# ---------------------------------------------------------------------------
# UniProt functional description (protein name, function, keywords, review status)
# ---------------------------------------------------------------------------

def _parse_protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    if rec:
        return str(rec.get("fullName", {}).get("value", "")).strip()
    subs = desc.get("submissionNames", [])
    if subs:
        return str(subs[0].get("fullName", {}).get("value", "")).strip()
    return ""


def _uniprot_functional_get(params: dict, verbose: bool = False) -> "list | None":
    if genomics_ora.HTTP_SESSION is None:
        return None
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            resp = genomics_ora.HTTP_SESSION.get(
                UNIPROT_SEARCH_URL, params=params, headers=genomics_ora.OMISCERA_HEADERS, timeout=DEFAULT_TIMEOUT,
            )
            if resp.status_code == 200:
                return resp.json().get("results", [])
            last_error = f"HTTP {resp.status_code}"
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < DEFAULT_RETRIES:
                time.sleep(1.5)
                continue
            break
        except Exception as exc:
            last_error = str(exc)
            if attempt < DEFAULT_RETRIES:
                time.sleep(1.0)
    if verbose:
        print(f"[GeneFunctionAnnotation] UniProt request failed: {last_error}")
    return None


def uniprot_accession_exists(accession: str, verbose: bool = False) -> bool:
    """True if UniProt recognizes this accession. Used by single-gene mode
    before skipping resolution for something that merely LOOKS like a UniProt
    accession (matches the regex) - a typo (e.g. P99999) or, rarely, a gene
    symbol that happens to match the accession pattern would otherwise be
    treated as a real accession, fetch nothing, and silently show up as a
    blank "Sparse" row instead of a clear error or a correct symbol lookup."""
    entries = _uniprot_functional_get(
        {"query": f"accession:{accession}", "format": "json", "size": 1, "fields": "accession"},
        verbose=verbose,
    )
    return bool(entries)


def fetch_uniprot_functional(
    accessions: list,
    batch_size: int = UNIPROT_BATCH_SIZE,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Fetch UniProt's functional description, keywords, and review status
    per accession — fields the shared gene resolver does NOT request (it only
    fetches accession/name/gene/proteome for resolution), so this is a
    separate fetch and does not touch Genomics_module_retriever.UNIPROT_FIELDS
    (ORA/Burden/Netprop's resolution calls stay lean).

    Returns {accession: {"protein_name", "organism", "function_text",
    "keywords": [str], "reviewed": bool, "protein_existence": str}}.
    function_text is the FIRST FUNCTION comment's prose, taken from UniProt's
    structured JSON (comments[].texts[].value) — already free of the
    flat-file's inline {ECO:...} evidence tags, since those are carried
    separately in each text's own "evidences" list, not baked into the string.
    """
    result: dict = {acc: _UNIPROT_FUNCTIONAL_CACHE[acc] for acc in accessions if acc in _UNIPROT_FUNCTIONAL_CACHE}
    to_fetch = [acc for acc in dict.fromkeys(accessions) if acc not in _UNIPROT_FUNCTIONAL_CACHE]
    total = len(to_fetch)
    done = 0

    for i in range(0, len(to_fetch), batch_size):
        batch = to_fetch[i:i + batch_size]
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        params = {"query": query, "format": "json", "size": len(batch), "fields": UNIPROT_FUNCTIONAL_FIELDS}
        entries = _uniprot_functional_get(params, verbose=verbose) or []

        seen = set()
        for entry in entries:
            acc = str(entry.get("primaryAccession", "")).strip()
            if not acc:
                continue
            seen.add(acc)
            func_comments = [c for c in (entry.get("comments") or []) if c.get("commentType") == "FUNCTION"]
            function_text = ""
            if func_comments:
                texts = func_comments[0].get("texts") or []
                if texts:
                    function_text = str(texts[0].get("value", "")).strip()
            keywords = [str(k.get("name", "")) for k in (entry.get("keywords") or []) if k.get("name")]
            entry_type = str(entry.get("entryType", ""))
            info = {
                "protein_name": _parse_protein_name(entry),
                "organism": str(entry.get("organism", {}).get("scientificName", "")),
                "function_text": function_text,
                "keywords": keywords,
                "reviewed": entry_type.startswith("UniProtKB reviewed"),
                "protein_existence": str(entry.get("proteinExistence", "")),
            }
            result[acc] = info
            _UNIPROT_FUNCTIONAL_CACHE[acc] = info

        for acc in batch:
            if acc not in seen:
                empty = {
                    "protein_name": "", "organism": "", "function_text": "",
                    "keywords": [], "reviewed": False, "protein_existence": "",
                }
                result[acc] = empty
                _UNIPROT_FUNCTIONAL_CACHE[acc] = empty

        done += len(batch)
        if progress_callback:
            try:
                progress_callback(done, total, f"Fetching UniProt functional data ({done}/{total} genes)")
            except Exception:
                pass

    return {acc: result.get(acc, {}) for acc in accessions}


# ---------------------------------------------------------------------------
# Annotation-confidence classification (the "dark gene" indicator)
# ---------------------------------------------------------------------------

def classify_annotation_confidence(go_terms: list, uniprot_info: dict) -> dict:
    """Classify how well-characterized a gene's annotation is, from its FULL
    GO evidence (all evidence codes, unfiltered by any display toggle) and
    UniProt curation status. Three tiers, in order of evidence strength:

      Well-characterized: reviewed (Swiss-Prot) AND >=1 experimental GO
        annotation AND a curated UniProt function description exist.
      Moderate: reviewed with SOME non-electronic support (an experimental GO
        term, a curated non-experimental GO term (e.g. ISS/TAS/IC), or a
        function description) but not meeting the "well-characterized" bar;
        OR unreviewed (TrEMBL) but with at least one experimental GO term.
      Sparse ("dark gene"): everything else — annotations are electronic
        (IEA) only and/or there is no function description and/or the entry
        is unreviewed with no experimental support. These are the genes most
        likely to reward further literature/experimental follow-up, since
        current databases have the least direct evidence for them.

    Returns {"status", "label", "n_experimental_terms", "n_total_terms",
    "reviewed", "has_function"} — the signals are returned alongside the
    verdict so the classification is always explainable, not a black box.
    """
    n_experimental = sum(1 for t in go_terms if t.get("evidence") in EXPERIMENTAL_EVIDENCE_CODES)
    n_curated_nonexp = sum(1 for t in go_terms if t.get("evidence") in CURATED_NONEXPERIMENTAL_EVIDENCE_CODES)
    n_total = len(go_terms)
    reviewed = bool(uniprot_info.get("reviewed"))
    has_function = bool(uniprot_info.get("function_text"))

    if reviewed and n_experimental >= 1 and has_function:
        status = "well_characterized"
    elif (reviewed and (n_experimental >= 1 or n_curated_nonexp >= 1 or has_function)) or (not reviewed and n_experimental >= 1):
        status = "moderate"
    else:
        status = "sparse"

    return {
        "status": status,
        "label": STATUS_LABELS[status],
        "n_experimental_terms": n_experimental,
        "n_total_terms": n_total,
        "reviewed": reviewed,
        "has_function": has_function,
    }


# ---------------------------------------------------------------------------
# Orchestration: fetch + classify + tabulate
# ---------------------------------------------------------------------------

def run_gene_function_annotation(
    gene_accession_pairs: list,
    verbose: bool = False,
    progress_callback=None,
) -> dict:
    """Full per-gene annotation pipeline. gene_accession_pairs: [(input_gene_id,
    accession), ...] (accessions may repeat across genes only in odd edge
    cases; duplicates are fetched once).

    ALWAYS returns full evidence (including IEA) in both tables — there is no
    include_iea parameter here on purpose: the "Include IEA" toggle on the
    results screen is a pure DISPLAY filter applied client-side (by grouping
    the tidy table), not a re-fetch. That is what makes it "live" (instant,
    no server round-trip) rather than requiring a fresh run for every toggle
    flip. Filtering evidence at this layer would throw away exactly the rows
    the toggle needs to reconstruct the "IEA excluded" view.

    Returns {"collapsed": DataFrame (one row/gene, full-evidence term lists),
    "tidy": DataFrame (one row/gene-term, ALL evidence), "counts":
    {well_characterized, moderate, sparse, n_genes}}.
    """
    accessions = [acc for _, acc in gene_accession_pairs if acc]
    unique_accessions = list(dict.fromkeys(accessions))

    def _go_progress(current, total, label):
        if progress_callback:
            progress_callback(current, total * 2, label)

    def _uniprot_progress(current, total, label):
        if progress_callback:
            progress_callback(total + current, total * 2, label)

    go_by_acc = fetch_go_annotations(unique_accessions, progress_callback=_go_progress, verbose=verbose)
    uniprot_by_acc = fetch_uniprot_functional(unique_accessions, progress_callback=_uniprot_progress, verbose=verbose)

    collapsed_rows, tidy_rows = [], []
    counts = {"well_characterized": 0, "moderate": 0, "sparse": 0}

    for gene_id, acc in gene_accession_pairs:
        terms = go_by_acc.get(acc, [])
        info = uniprot_by_acc.get(acc, {})
        classification = classify_annotation_confidence(terms, info)
        counts[classification["status"]] += 1

        by_ns: dict = {"biological_process": [], "molecular_function": [], "cellular_component": []}
        for t in terms:
            ns = t.get("namespace")
            if ns in by_ns and t.get("name"):
                by_ns[ns].append(t["name"])
            tidy_rows.append({
                "Gene": gene_id, "UniProt ID": acc, "GO ID": t.get("go_id", ""),
                "Term": t.get("name", ""), "Namespace": ns or "",
                "Evidence Code": t.get("evidence", ""), "Qualifier": t.get("qualifier", ""),
            })

        collapsed_rows.append({
            "Gene": gene_id,
            "UniProt ID": acc,
            "Protein Name": info.get("protein_name", ""),
            "Organism": info.get("organism", ""),
            "Reviewed": bool(info.get("reviewed")),
            "GO Biological Process": "; ".join(sorted(set(by_ns["biological_process"]))),
            "GO Molecular Function": "; ".join(sorted(set(by_ns["molecular_function"]))),
            "GO Cellular Component": "; ".join(sorted(set(by_ns["cellular_component"]))),
            "Keywords": "; ".join(info.get("keywords") or []),
            "UniProt Function": info.get("function_text", ""),
            "Annotation Status": classification["label"],
            "Experimental GO Terms": classification["n_experimental_terms"],
            "Total GO Terms": classification["n_total_terms"],
        })

    collapsed_df = pd.DataFrame(collapsed_rows, columns=ANNOTATION_COLLAPSED_COLUMNS)
    tidy_df = pd.DataFrame(tidy_rows, columns=ANNOTATION_TIDY_COLUMNS)
    counts["n_genes"] = len(gene_accession_pairs)
    return {"collapsed": collapsed_df, "tidy": tidy_df, "counts": counts}
