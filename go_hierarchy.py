"""
Shared GO (Gene Ontology) OBO hierarchy parsing, used by both
Genomics_module.py and Transcriptomics_module.py so a gene annotated to a
specific term also counts toward that term's ancestors (is_a / part_of),
matching clusterProfiler/goseq/topGO/DAVID - none of which restrict
enrichment to only the most specific (leaf-level) annotation a gene happens
to carry. QuickGO's downloadSearch/annotation APIs return DIRECT annotations
only, so without this, every non-leaf GO term is systematically
under-populated relative to those reference tools.

Both modules download/cache the identical go-basic.obo file (same URL, same
CACHE_DIR), so this module parses it once and both reuse the result.
"""
import os

# path -> (mtime, parsed_dict); avoids re-parsing the ~2M-line OBO file on
# every ORA run within the same server process.
_PARSED_CACHE: dict = {}


def parse_go_obo(path: str, verbose: bool = False) -> dict:
    """Parse a go-basic.obo file into term metadata plus its ancestor closure.

    Returns:
        {
          "names": {go_id: name},              # includes alt_id aliases
          "namespace": {go_id: namespace_str},  # e.g. "biological_process"
          "ancestors": {go_id: frozenset(ancestor_go_id, ...)},  # transitive
              is_a/part_of closure, restricted to same-namespace edges only
              (go-basic.obo's own edges are already same-namespace by
              construction; this is a defensive guard, not an assumption).
          "alt_ids": {alt_go_id: canonical_go_id},
        }

    Obsolete terms (is_obsolete: true) are dropped entirely - they never
    appear as keys in any of the four dicts, so they cannot appear in
    results or receive/contribute propagated genes.
    """
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        mtime = None
    cached = _PARSED_CACHE.get(path)
    if cached is not None and cached[0] == mtime:
        return cached[1]

    names: dict = {}
    namespace_of: dict = {}
    raw_parents: dict = {}
    alt_ids: dict = {}

    cur_id, cur_name, cur_ns = "", "", ""
    cur_alts: list = []
    cur_parents: set = set()
    cur_obsolete = False
    in_term = False

    def _flush():
        if not cur_id or cur_obsolete:
            return
        names[cur_id] = cur_name or cur_id
        namespace_of[cur_id] = cur_ns
        raw_parents[cur_id] = set(cur_parents)
        for a in cur_alts:
            alt_ids[a] = cur_id
            names.setdefault(a, names[cur_id])  # alias resolves to the canonical term's name

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            # go-basic.obo's trailing [Typedef] stanzas (has_ontology_root_term,
            # part_of, regulates, ...) also carry `name:`/`namespace:` lines.
            # Flushing/resetting on ANY "[...]" header - not just "[Term]" - and
            # only accumulating fields while inside a [Term] stanza prevents
            # those from overwriting the last real term's state (verified: the
            # last GO term in the file was previously corrupted to name
            # "regulates", namespace "external", with its ancestors dropped).
            if line.startswith("[") and line.endswith("]"):
                _flush()
                cur_id, cur_name, cur_ns = "", "", ""
                cur_alts, cur_parents = [], set()
                cur_obsolete = False
                in_term = (line == "[Term]")
                continue
            if not in_term:
                continue
            if line.startswith("id: GO:"):
                cur_id = line[4:].strip()
            elif line.startswith("name:"):
                cur_name = line[5:].strip()
            elif line.startswith("namespace:"):
                cur_ns = line[10:].strip()
            elif line.startswith("alt_id: GO:"):
                cur_alts.append(line[8:].strip())
            elif line.startswith("is_a: GO:"):
                # "is_a: GO:0008150 ! biological_process" - keep just the ID.
                cur_parents.add(line[6:].strip().split()[0])
            elif line.startswith("relationship: part_of GO:"):
                cur_parents.add(line[len("relationship: part_of "):].strip().split()[0])
            elif line.startswith("is_obsolete: true"):
                cur_obsolete = True
        _flush()

    # Resolve any parent referenced by an alt_id (not yet re-pointed to its
    # canonical id) and drop dangling/obsolete parent references.
    def _canon(go_id: str) -> str:
        return alt_ids.get(go_id, go_id)

    parents = {
        gid: {_canon(p) for p in ps if _canon(p) in names}
        for gid, ps in raw_parents.items()
    }

    ancestors: dict = {}

    def _closure(gid: str, seen: frozenset) -> frozenset:
        if gid in ancestors:
            return ancestors[gid]
        if gid in seen:
            return frozenset()  # defensive cycle guard; GO is a DAG
        seen = seen | {gid}
        acc: set = set()
        for p in parents.get(gid, ()):
            if namespace_of.get(p) != namespace_of.get(gid):
                continue
            acc.add(p)
            acc |= _closure(p, seen)
        result = frozenset(acc)
        ancestors[gid] = result
        return result

    for gid in parents:
        _closure(gid, frozenset())

    parsed = {"names": names, "namespace": namespace_of, "ancestors": ancestors, "alt_ids": alt_ids}
    if verbose:
        print(f"[GO hierarchy] parsed {len(names)} terms, {sum(1 for a in ancestors.values() if a)} with ancestors")
    _PARSED_CACHE[path] = (mtime, parsed)
    return parsed


def propagate_term_map(term_map: dict, hierarchy: dict) -> dict:
    """Expand a {go_id: {"name": str, "genes": set|list}} map of DIRECT
    annotations (one namespace's worth) so each gene also counts toward every
    ancestor of its directly-annotated term, transitively.

    Terms that only exist as ancestors (never directly annotated in
    term_map) are added as new entries, named from the hierarchy. Safe to
    call with terms not present in `hierarchy` (e.g. a stale/obsolete id from
    a live QuickGO response racing an OBO refresh) - they pass through with
    no propagation rather than raising.

    An incoming id that is a secondary/alt_id is canonicalized first, so its
    genes merge into the canonical term (and inherit its ancestors) instead
    of sitting as an orphan duplicate with no propagation.
    """
    ancestors = hierarchy["ancestors"]
    names = hierarchy["names"]
    alt_ids = hierarchy.get("alt_ids", {})
    expanded: dict = {}
    for raw_go_id, info in term_map.items():
        go_id = alt_ids.get(raw_go_id, raw_go_id)
        genes = set(info["genes"])
        bucket = expanded.setdefault(go_id, {"name": info.get("name") or names.get(go_id, go_id), "genes": set()})
        bucket["genes"] |= genes
        for anc in ancestors.get(go_id, ()):
            anc_bucket = expanded.setdefault(anc, {"name": names.get(anc, anc), "genes": set()})
            anc_bucket["genes"] |= genes
    return expanded
