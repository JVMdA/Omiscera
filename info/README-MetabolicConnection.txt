METABOLIC CONNECTION
=======================

WHAT THIS MODULE IS FOR
--------------------------
You give it a list of compound names (or, optionally, identifiers you already
have — PubChem ID, SMILES, InChIKey, HMDB ID, KEGG ID). For each compound, it
asks two separate questions:

  1. RESOLUTION — what is this compound? It looks the name up online and
     fills in every identifier it can find: PubChem ID, SMILES, InChIKey,
     HMDB ID, Formula, and KEGG ID.

  2. EVIDENCE — is there evidence this compound occurs in, or is metabolised
     by, a specific organism? You pick one or more databases and, for each,
     a biological context (a taxon, a species, or — for KEGG — an organism
     code), and the module checks each compound against it.

These are two genuinely separate steps, run in that order, because a
compound has to be identified before it can be looked up anywhere.

STAGE 1 — COMPOUND RESOLUTION
----------------------------------
Upload a table with a "Name" column (optionally alongside any identifiers
you already have — a header named "KEGG ID"/"KEGGID"/"KEGG" is recognised
automatically, same as the other identifier columns).

Each compound's Mapping_status is:
  - mapped   — PubChem ID, SMILES, InChIKey, HMDB ID AND KEGG ID were all
               resolved.
  - partial  — at least the name resolved to something, but one or more of
               those five identifiers is still missing.
  - not_found — nothing could be resolved.

Both "mapped" and "partial" compounds proceed to Stage 2 — the label only
changes what is DISPLAYED and counted, never what gets analysed. A compound
missing only its KEGG ID, for example, still gets a full Rhea + UniProt /
LOTUS / HMDB / NPAtlas evidence check; it just cannot appear in the KEGG
enzyme table (Table 3), since that table is keyed on KEGG ID.

Resolved results are cached (by name + any identifiers you supplied), so
re-running the same compound list does not re-hit any external API.

STAGE 2 — EVIDENCE CHECK, PER DATABASE
-------------------------------------------
Five databases are available, each checked independently per rank/context you
configure:

  - LOTUS      — natural-product occurrence records, searched by taxon.
  - HMDB       — human metabolome evidence, either searched online or from an
                 HMDB export you upload.
  - Rhea + UniProt — is the compound a substrate/product of a reaction whose
                 enzyme exists (reviewed or unreviewed) in the chosen UniProt
                 taxon? Produces Table 2 (enzyme presence).
  - NPAtlas    — microbial (bacteria/fungi) natural-product occurrence,
                 searched by genus/family/species.
  - KEGG       — see below. Produces Table 3 (enzyme presence), a separate
                 table from Rhea + UniProt's because KEGG's data supports an
                 extra column (Pathway) that Rhea's does not.

LOTUS, HMDB, Rhea + UniProt and NPAtlas results land in Table 1 (a simple
Found/Not_found/Not_evaluated grid, one column per database) except Rhea +
UniProt, which — like KEGG — is "enzyme-shaped" and gets its own table
instead of a single grid cell, because "does this compound have evidence in
this organism" for an enzyme-linked database really means "which enzymes,
and is each one present", not a single yes/no.

KEGG, SPECIFICALLY
-----------------------
KEGG answers the same enzyme-presence question as Rhea + UniProt, but through
KEGG's own data and REST API, which is why it produces its own table (Table
3) rather than being merged into Table 2.

How a KEGG row is built, for one compound in one organism:
  1. The compound's KEGG ID's linked enzymes come from a DIRECT KEGG link
     (compound -> enzyme), not by first finding reactions and then looking at
     their enzymes.
  2. Substrate/product ROLE is a SEPARATE lookup: the compound's reactions are
     found, each reaction's equation is parsed on its "<=>" split to see
     which side the compound is on, and each reaction's enzymes are pulled in
     turn. This is joined back to the enzyme list from step 1 by EC number.
     An enzyme from step 1 with no role found this way is labelled "Unknown"
     — it is never dropped just because the role lookup didn't cover it.
  3. PRESENCE is simply: is this EC number in the chosen organism's KEGG gene
     set (one lookup per organism, reused across every compound in that
     rank, not refetched per compound)?
  4. PATHWAY lists the named KEGG pathways (excluding KEGG's generic
     "overview" maps, which are not real organism-specific pathways) that
     enzyme participates in, within that organism.

KEGG ID RESOLUTION. A compound's KEGG ID is found by searching KEGG by name,
which can return more than one candidate for an ambiguous name. Each
candidate is VERIFIED, not just name-matched: KEGG's own record cross-
references PubChem, but that cross-reference is a PubChem Substance ID (SID),
not a Compound ID (CID) — a legacy KEGG/NIH data-submission detail. The SID is
resolved to its real CID via PubChem's own substance-to-compound endpoint,
and a candidate is only accepted if that CID matches the PubChem ID this
module already resolved independently for the same compound. This matters:
treating the SID as if it were a CID directly would silently match a
different, unrelated molecule for some compounds, since the two ID spaces
overlap numerically without meaning the same thing. Only if no candidate
carries a PubChem cross-reference at all does the module fall back to the
best unverified name match.

REQUEST BEHAVIOUR. Batched lookups (used for the Role column's reaction
equations) request up to 10 entries per call — KEGG's own documented maximum
— rather than one request per reaction. All KEGG requests are throttled to
roughly 3 per second, matching KEGG's published usage guidance. Every KEGG
response is parsed by matching its own ENTRY line back to the ID it belongs
to, not by request position — KEGG's batched endpoint silently DROPS invalid
IDs from its response rather than returning a placeholder for them, so
matching by position would risk attributing one entry's data to a different
ID whenever a batch contained even one bad ID.

WHY KEGG, GIVEN THE OTHER FOUR DATABASES. For non-model organisms in
particular, KEGG's enzyme/pathway coverage is materially more complete than
the free alternatives this module also offers. It was removed from an
earlier version of this app and is being re-added now.

ACADEMIC / LICENSING NOTE. KEGG's REST API is intended for individual,
academic use, and its terms restrict systematic bulk downloading and
commercial redistribution of its data. This module's queries are targeted,
per-compound, interactive lookups driven by what you actually upload — the
same pattern KEGG's REST API exists to support — not a bulk export of KEGG's
database. This app is distributed as source for you to run yourself, rather
than offered as a hosted service making KEGG calls on other users' behalf.

RESULTS AND DOWNLOAD
-------------------------
Per rank, you get:
  - Table 1 — Found/Not_found/Not_evaluated per non-enzyme database.
  - Table 2 — Rhea + UniProt enzyme presence (Name, Enzyme, Role, Enzyme
    presence), sorted Found first.
  - Table 3 — KEGG enzyme presence (Name, Enzyme, Role, Enzyme presence,
    Pathway), sorted Found first.

Each table downloads independently as CSV.
