"""Live smoke test: hit the real TRRUST/OmniPath/ChEA3 endpoints
through Transcriptomics_regulon's actual fetchers (no mocking) and sanity-
check the parsed shapes against well-known TF biology."""
import os
import sys

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_regulon as tr

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


print("=== TRRUST (human) live fetch ===")
trrust = tr.fetch_trrust_regulons("9606", force_refresh=True, verbose=True)
check(not trrust["fetch_failed"], f"fetch did not fail (note={trrust.get('note')})")
check(len(trrust["regulons"]) > 500, f"parsed a plausible number of TFs (got {len(trrust['regulons'])})")
tp53 = trrust["regulons"].get("TP53")
check(tp53 is not None, "TP53 present (a famous, heavily-studied TF)")
if tp53:
    check(len(tp53["genes"]) > 50, f"TP53 has a large target list (got {len(tp53['genes'])})")
    check(len(tp53["weights"]) > 0, f"TP53 has at least some signed edges (got {len(tp53['weights'])})")

print("=== DoRothEA (human, confidence A) live fetch ===")
doro = tr.fetch_dorothea_regulons("9606", confidence_levels="A", force_refresh=True, verbose=True)
check(not doro["fetch_failed"], f"fetch did not fail (note={doro.get('note')})")
check(len(doro["regulons"]) > 20, f"parsed a plausible number of high-confidence TFs (got {len(doro['regulons'])})")
myc = doro["regulons"].get("MYC")
check(myc is not None, "MYC present at confidence A (a canonical high-confidence TF)")

print("=== DoRothEA (mouse) live fetch ===")
doro_mouse = tr.fetch_dorothea_regulons("10090", confidence_levels="AB", force_refresh=True, verbose=True)
check(not doro_mouse["fetch_failed"], f"mouse fetch did not fail (note={doro_mouse.get('note')})")
check(len(doro_mouse["regulons"]) > 5, f"parsed some mouse TFs (got {len(doro_mouse['regulons'])})")

print("=== CollecTRI (human) live fetch - complex exclusion ===")
collectri = tr.fetch_collectri_regulons("9606", force_refresh=True, verbose=True)
check(not collectri["fetch_failed"], f"fetch did not fail (note={collectri.get('note')})")
complex_like = [tf for tf in collectri["regulons"] if "_" in tf]
check(len(complex_like) == 0, f"no COMPLEX-derived multi-gene TF names leaked through (got {complex_like[:5]})")
jun = collectri["regulons"].get("JUN")
check(jun is not None, "JUN present in CollecTRI")

print("=== upstream-coverage watch: TRRUST still publishes only human+mouse ===")
# TRRUST_URL_BY_TAXON (Transcriptomics_regulon.py) hardcodes exactly two
# organism files because that is genuinely all TRRUST has ever published.
# This is the mechanism that makes "if TRRUST adds another organism, will
# the app catch it?" actually true: probe realistic next-candidate organism
# names against TRRUST's own real URL pattern and confirm they still 404.
# The moment one of these stops 404ing, it's time to widen
# TRRUST_URL_BY_TAXON, not something that should go unnoticed.
_TRRUST_CANDIDATE_ORGANISMS = ["rat", "zebrafish", "fly", "drosophila", "worm", "yeast"]
if tr.tx.HTTP_SESSION is not None:
    still_404 = []
    unexpected = []
    for org in _TRRUST_CANDIDATE_ORGANISMS:
        url = f"https://www.grnpedia.org/trrust/data/trrust_rawdata.{org}.tsv"
        try:
            resp = tr.tx.HTTP_SESSION.get(url, headers=tr.tx.OMISCERA_HEADERS, timeout=30)
            if resp.status_code == 404:
                still_404.append(org)
            else:
                unexpected.append((org, resp.status_code))
        except Exception as exc:
            unexpected.append((org, str(exc)))
    check(len(unexpected) == 0,
          f"every candidate-organism URL TRRUST doesn't (yet) support still 404s as expected - a non-404 "
          f"here means TRRUST just published a new organism file (got unexpected: {unexpected})")
    check(set(still_404) == set(_TRRUST_CANDIDATE_ORGANISMS),
          f"sanity: all {len(_TRRUST_CANDIDATE_ORGANISMS)} candidates were actually reachable and checked "
          f"(got {len(still_404)} confirmed 404s)")
    # SENSITIVITY: prove the detection logic really does flag a non-404 as
    # "unexpected", using a REAL live response rather than a fabricated one -
    # probe "human", which genuinely exists (200), through the exact same
    # candidate-checking logic above.
    sentinel_unexpected = []
    resp = tr.tx.HTTP_SESSION.get(
        "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv", headers=tr.tx.OMISCERA_HEADERS, timeout=30,
    )
    if resp.status_code != 404:
        sentinel_unexpected.append(("human", resp.status_code))
    check(len(sentinel_unexpected) == 1 and sentinel_unexpected[0][1] == 200,
          f"the SAME detection logic correctly flags a real, live, genuinely-existing organism file "
          f"('human', which TRRUST does publish) as non-404 - proving the watch above would truly catch "
          f"it if a candidate organism stopped 404ing (got {sentinel_unexpected})")
else:
    check(False, "requests/HTTP_SESSION not available - could not run the upstream-coverage watch")

print("=== ChEA3 (one library, to keep this fast) live fetch ===")
chea3 = tr.fetch_chea3_regulons(libraries=["ENCODE_ChIP-seq"], force_refresh=True, verbose=True)
check(not chea3["fetch_failed"], f"fetch did not fail (note={chea3.get('note')})")
check(chea3["signed"] is False, "ChEA3 flagged unsigned")
stat1 = chea3["regulons"].get("STAT1")
check(stat1 is not None, "STAT1 present in ENCODE_ChIP-seq")
if stat1:
    check(len(stat1["genes"]) > 100, f"STAT1 has a large merged target list across cell lines (got {len(stat1['genes'])})")
    check(stat1["weights"] == {}, "STAT1 weights empty (unsigned)")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
