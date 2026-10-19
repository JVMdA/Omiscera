===========================================================================
OMISCERA TEST SUITE - Transcriptomics, Genomics and Proteomics
===========================================================================

These scripts were written while building and auditing the Transcription
Factor / Regulon module, then extended when MSigDB was added to ORA/GSEA,
again for Transcriptomics Network / Module Detection, and again for the
Proteomics module (Over-Representation Analysis). A few also cover Genomics
screens and shared, app-wide code (the table renderer, organism-state
guards, the Python<->JS allowlist mirrors).
They are plain Python scripts, NOT pytest: each one
does its own PASS/FAIL bookkeeping, prints one line per check, and exits 0
(all passed) or 1 (something failed). Run any of them directly:

    python tests\test_regulon_engine.py

Every script is repo-relative - it locates the modules under test from its
own file location - so the suite runs from any checkout without editing
paths.

WHY THEY EXIST IN THIS FORM. The module's statistics could not be verified
by inspection alone, so essentially every claim it makes is checked
numerically here rather than asserted: the ULM closed form against a
brute-force OLS fit, the sign conventions end to end, the empirical-null
estimator against simulated ground truth, and each figure's encoding against
the data it claims to show. Several of these scripts also contain
SENSITIVITY checks that deliberately break the implementation and confirm
the corresponding test goes red - a test that cannot fail is not evidence.


---------------------------------------------------------------------------
NO SERVER NEEDED (pure functions, hand-built fixtures)
---------------------------------------------------------------------------
test_regulon_engine.py         117 checks. The core suite: ULM closed form,
                               sign correctness, degenerate inputs, parsers,
                               cache round-trip, and the empirical-null
                               calibration (section 12b covers the
                               truncated-normal estimator, including that it
                               still recovers GENUINE inflation - the
                               condition that matters most - plus a
                               sensitivity check that breaks the truncation
                               correction and confirms failure).

test_example_files.py           22 checks. The bundled example DE files parse,
                                resolve, and carry the columns the app expects.

audit_figures_and_math.py       16 checks. Full figure + mathematics audit:
                                cross-figure consistency, coherence-bar signs,
                                QQ-plot quantiles, empirical-null behaviour as
                                the truly-active fraction rises (including its
                                documented breakdown past ~55%), p-value
                                underflow, and multiple-testing families.

audit_lambda_ambiguity.py        4 checks. Where the lambda ambiguity now
                                lives after the truncated-normal estimator
                                replaced plain MAD. Records a DISPROVEN
                                hypothesis (the QQ plot's central slope does
                                not separate inflation from signal) rather
                                than dropping it.

audit_regulon_untested.py       15 checks. Sweeps for module surface that no
                                other test touches.

verify_regression_fix.py         5 checks. Guards the fix to the ULM
                                regression figure (it once plotted log2FC
                                against a slope fitted in ranking-metric
                                units).

smoke_calibration_wiring.py     16 checks. The calibrated columns the engine
                                emits really are recomputable from the stored
                                null centre and lambda.

probe_better_null_estimator.py  Feasibility evidence, kept for provenance:
probe_truncated_validity.py     the measurements that justified replacing MAD
calibration_regulon_activity.py with the truncated-normal MLE, and the
                                original real-data measurement that justified
                                having calibration at all. Cited by name in
                                Transcriptomics_regulon.py's docstrings.
                                These PRINT tables rather than asserting; the
                                assertions were promoted into
                                test_regulon_engine.py section 12b.

calibration_network_production_scale.py  Feasibility/calibration evidence
                                for the network engine's skew-aware p-value
                                and permutation-count defaults, same
                                provenance-script precedent as the group
                                above (prints tables, no PASS/FAIL
                                bookkeeping, real runtime ~35 min - not part
                                of the routine suite). Part A: true-null
                                calibration at PRODUCTION scale (16,201
                                proteins, the real cached human STRING graph,
                                not the 60-node random graphs
                                test_network_engine_core.py's own calibration
                                uses) - 5 repeats of random, structure-free
                                seeds, so every FDR<0.05 call is by
                                construction a false positive. Measured: the
                                skew-aware p-value calls ZERO false positives
                                across all 5 repeats; the pre-fix normal-tail
                                approximation, recomputed from the same null
                                moments (no second permutation run needed),
                                calls ~3.4% (~551/16,201) significant in
                                EVERY repeat - several hundred false positives
                                on a genuinely null input, every single run.
                                Part B: the 200/500/1000/2000 permutation-
                                count sweep on the bundled human example that
                                justified dropping 200 and raising the
                                default to 1000 - see both network READMEs
                                and the Fix 3 commit message for the full
                                finding (a specific protein's null_std jumped
                                10x between 200 and 500 permutations, a real
                                null-variance-underestimation artifact, not
                                imprecision).

generate_full_examples.py       Regenerates the bundled Data3 example DE
                                files. Not a test - run only to rebuild them.

test_table_long_list_rendering.py  10 checks. A real, user-reported bug in
                                buildTable() (static/app.js) - the SHARED
                                table renderer used by EVERY results screen
                                in the app, not MSigDB-specific code. A long
                                comma-joined gene-list cell (Overlapping
                                Genes/Leading Edge Genes/Genes) wrapping in a
                                narrow column inflated its whole ROW to
                                thousands of pixels (3295px, measured live);
                                since cells default to vertical-align:middle,
                                every OTHER cell in that row - everything a
                                user actually scans - was centered far
                                outside the results table's scroll viewport,
                                making a table with real data look completely
                                empty. Fixed with a 3-line CSS clamp on those
                                specific columns (full list still in the
                                title tooltip and the CSV download) plus
                                vertical-align:top as an independent
                                safeguard. A sensitivity check proves the
                                row-height assertion is not vacuous by
                                measuring the same input with the clamp
                                withheld.

test_organism_consistency.py    24 checks. Two pre-existing organism-state
                                bugs, older and more dangerous than anything
                                MSigDB introduced, found by a post-
                                implementation audit: taxon_name used to
                                commit synchronously in the resolve ENDPOINT
                                while taxon_id only committed inside the
                                resolve JOB on success, so a failed re-resolve
                                left the two disagreeing (taxon_id drives GO/
                                MSigDB/gene-length, taxon_name drives Reactome
                                - a drifted session doesn't error, it silently
                                mixes two organisms in one result set). Tests,
                                with hand-built session dicts and direct
                                in-process calls into API.py (no server, no
                                network): _require_consistent_taxon's own
                                logic; that _job_resolve_de now commits both
                                fields together, atomically, with a
                                sensitivity check proving a failed resolve
                                under the OLD ordering really did leave them
                                disagreeing; and that the guard is actually
                                WIRED into every transcriptomics configure
                                endpoint - ORA/GSEA/regulon/network, DISCOVERED
                                from API.py's own "@app.post(.../configure)"
                                decorators rather than hand-listed, so a future
                                5th module's route is covered automatically
                                the day it exists. The network route's
                                consistent-session check stubs
                                genomics_net.is_taxon_in_string so this file
                                keeps its "no live server, no network" promise
                                - the drifted-session branch never reaches
                                that call for any route, since the guard fires
                                first.

test_network_engine_core.py     58 checks. Characterization tests for the
                                SHARED numerical engine in Genomics_network.py
                                (STRING parsing, random-walk-with-restart, the
                                degree-preserving permutation null, MCL) -
                                until this file, ZERO unit coverage existed
                                for any of it, discovered while building the
                                Transcriptomics Network / Module Detection
                                module directly on top of it. Hand-built
                                gzipped STRING fixture files (no real
                                downloads): the confidence-cutoff boundary
                                (>=, not >, proven with two real calls at 700
                                vs 701), propagate_rwr against a closed-form
                                sparse solve, the permutation null's
                                determinism AND degree preservation (a spy
                                intercepts every permuted seed vector and
                                compares its degree-bin composition against
                                the real seeds'), a CHARACTERIZED limitation
                                (a graph under ~20 proteins makes every degree
                                bin hold exactly one candidate, so the null
                                can never move and nothing can ever be
                                significant - asserted so this can't silently
                                change), FDR calibration under a TRUE null (5
                                repeats, random graphs, no planted signal,
                                worst call rate 1.7%), POWER and SPECIFICITY
                                together on a planted two-cluster fixture
                                (seeded cluster's satellites significantly
                                out-score an unseeded, degree-MATCHED control
                                cluster's satellites), mcl_cluster/
                                detect_modules, and a GOLDEN hit table the
                                Transcriptomics engine's extraction had to
                                reproduce byte-identically. Also covers three
                                evaluation-driven fixes: _GRAPH_CACHE is
                                thread-safe under a 40-thread hammering of
                                build_string_network across 5 keys (0 errors,
                                cache never exceeds its max size), with a
                                sensitivity check proving the new RLock is
                                actually held during mutation, not just
                                assumed thread-safe by luck; _welford_update's
                                extension to a third central moment (Pébay's
                                streaming skewness) matches a batch numpy
                                computation over genuinely right-skewed gamma
                                samples, and the resulting Pearson III p-value
                                tracks a planted fixture's own empirical
                                exceedance rate materially better than the
                                pre-fix normal-only p-value did, for the
                                highest-skew protein, with a sensitivity check
                                forcing pure-normal fallback (SKEW_GAMMA_
                                THRESHOLD=1e9) and confirming exactly the
                                Fix-1-specific checks go red; and a second,
                                independent recruited-only BH-FDR family
                                (fdrs_recruited_family) is NaN at every seed
                                position, matches an independent recomputation
                                to 1e-12, and genuinely differs from the
                                all-protein family's FDR for the same
                                proteins - confirmed to NOT match what BH over
                                all p-values (the pre-fix approach) would give.

test_transcriptomics_network.py 53 checks. Unit tests for
                                Transcriptomics_network.py, which had ZERO
                                unit coverage of its own (the shared engine
                                above is not retested here). The single
                                highest-value check: enrich_modules_by_symbol
                                uses the NETWORK's own gene symbols as its
                                enrichment universe, not the DE table's tested
                                genes - asserts a RECRUITED gene (absent from
                                the DE table by definition) is counted in k
                                and listed in "Overlapping Genes", with a
                                sensitivity check re-running the identical
                                real engine call against the DE-table-only
                                universe and watching the recruited gene
                                silently vanish (k drops 3 -> 2) - then
                                confirming the production code itself goes
                                red the same way when temporarily swapped to
                                match, and restored. Also: build_string_id_
                                index's collision resolution (best-connected
                                protein wins) against a naive dict(zip(...))
                                last-wins run on the SAME graph;
                                map_de_table_to_string's tiering (a true
                                partition, symbol beats a coincidentally-
                                matching Ensembl id) and its <10%
                                refusal / 10-50% warning gates;
                                build_seed_weights' equal weighting (two DEG
                                rows on one protein now genuinely weigh
                                exactly 1.0, not the old summed 2.0 - a real
                                fix, with a sensitivity check reverting the
                                production code to the old summing form and
                                confirming the check goes red), abs_metric's
                                sign-symmetry with a sensitivity check using
                                the real shift-positive branch to show what
                                breaks without it, its winsorization ratio
                                (exactly 20, by construction), and its own fix
                                (a multi-row-mapped protein takes the MAX
                                magnitude among its rows, not the sum,
                                discriminated against a hand-computed
                                sum-hypothesis value on real data); the
                                giant-module warning; the small-graph guard;
                                and enrich_modules_by_symbol's namespace-
                                selection correction (a two-namespace fixture
                                confirms top_fdr_corrected == min(raw_fdr *
                                n_namespaces_compared, 1.0) exactly, with two
                                sensitivity checks - a stand-in uncorrected
                                wrapper, and a real production-code revert -
                                both confirmed to go red). PLUS (added
                                alongside multi-database module labelling):
                                enrich_modules_by_symbol's `database` param
                                became `databases` (a list) - selecting GO +
                                MSigDB at once now returns BOTH candidates
                                under one module (`{"best", "candidates",
                                "n_namespaces_compared"}`, not just the
                                winner discarding the loser), with a fixture
                                using a FILLER term per database (the same
                                degenerate-universe trap FAKE_GO_SETS_TWO_NS
                                already worked around - a lone real term IS
                                the annotated universe, forcing p=1.0
                                regardless of overlap) so MSigDB's tighter
                                match genuinely outscores GO's looser one and
                                deterministically wins "best" while GO's own
                                candidate table stays independently usable.

test_transcriptomics_client_state.py  49 checks. Generic, source-level
                                invariants over static/app.js's Transcriptomics
                                wiring - same text-parsing technique as
                                test_organism_consistency.py, client-side. The
                                module list (ora/gsea/regulon/network) is
                                DISCOVERED from freshTranscriptomicsState()
                                itself, so a future 5th tile is covered here
                                automatically: every module's Config/Results/
                                RunLocked is genuinely cleared (plus its
                                _<mod>Draft deleted) on an organism change or
                                new upload, with a sensitivity check on an
                                in-memory copy missing one real line; every
                                module's router page routes to a function that
                                is actually DEFINED elsewhere in the file, not
                                just referenced; and every tile in the grid
                                array has either a not-yet-implemented flag or
                                a Continue-chain branch, with a sensitivity
                                check confirming a hypothetical unlisted tile
                                genuinely fails the same coverage check the 4
                                real tiles pass.

test_msigdb_fetcher.py          52 checks. The MSigDB fetcher added to
                                Transcriptomics_module.py: the availability
                                map (msigdb_collection_available - asserts
                                PID/Oncogenic have no mouse file, Cell type
                                does; a sensitivity check proves this isn't
                                vacuous by temporarily reintroducing the
                                c->m substitution bug it exists to prevent),
                                GMT parsing, the fetch_failed contract (HTTP
                                error/exception/unavailable-pair all degrade
                                safely without touching the cache), and that
                                the cache key includes MSIGDB_VERSION so a
                                version bump can't serve a stale entry. Fake
                                HTTP session, no real network.

test_reactome_module_organism_filter.py  8 checks. A real bug found while
                                auditing organism generality app-wide:
                                Genomics_network.enrich_modules (the
                                per-module functional-enrichment step behind
                                Network Propagation / Module Detection's
                                Reactome option) called
                                Genomics_module.build_reactome_universe
                                WITHOUT taxon_name. build_reactome_universe
                                filters by organism NAME (unlike GO, which
                                filters by taxon_id and is unaffected), so
                                the filter never ran - measured against the
                                real bulk file: 91,531 accessions pooled
                                across every organism vs 12,155 for human
                                alone, a >7x inflation that deflates fold
                                enrichment and corrupts p-values, worst for
                                non-model organisms. It was also cached under
                                the SAME (database, taxon_id) key every other
                                caller (Genomics ORA, Burden Test) uses for
                                that taxon, so one module-detection run could
                                poison Reactome results for every other
                                analysis on that taxon afterward. Fixed by
                                threading taxon_name through enrich_modules
                                and its API.py call site, plus a bumped
                                REACTOME_UNIVERSE_CACHE_VERSION so a pre-fix
                                poisoned cache entry is never re-served (a
                                fake poisoned entry is written under the OLD
                                key here and confirmed genuinely ignored, not
                                just assumed). Two sensitivity checks revert
                                the real fix (the taxon_name threading, then
                                the cache-version bump) and confirm each goes
                                red. Direct function calls against the real,
                                local Reactome bulk file - no server needed.

test_transcriptomics_kegg.py   78 checks. KEGG in Transcriptomics ORA +
                                GSEA, fully offline (every KEGG call
                                monkeypatched, on-disk cache redirected to a
                                temp dir so a run can neither read nor poison
                                the real cache - the mistake that had to be
                                fixed in test_kegg_gene_map.py). Transcriptomics
                                is gene-SYMBOL-native while KEGG is keyed by
                                KEGG gene ID, so KEGG needs two bridges nothing
                                else in that module needs, and each can fail
                                silently:
                                (1) ORGANISM. KEGG publishes no taxon ->
                                organism-code endpoint, so
                                resolve_kegg_organism_code matches on organism
                                NAME and then VERIFIES against KEGG's own
                                "TAXONOMY  TAX:<id>" line. Covers the exact
                                match, the genus+species fallback (the app's
                                "Saccharomyces cerevisiae S288C" vs KEGG's
                                unsuffixed genome name), and the case that
                                makes verification non-optional: a species-
                                level taxon whose only KEGG genomes are
                                strains, which must resolve to NOTHING rather
                                than silently test the wrong strain. This is
                                the same bug class as commit 27ce033 (KEGG
                                organism divergence in Genomics ORA), which is
                                why Transcriptomics derives the organism
                                instead of offering a second picker.
                                (2) IDENTIFIER. KEGG gene -> Entrez (KEGG's own
                                conv/ncbi-geneid table) -> this run's symbols
                                (the resolved DE table's own columns, now
                                shared with Reactome via
                                _entrez_sets_to_symbol_sets). Includes a
                                sensitivity check that a BLANK Entrez ID never
                                matches a blank entrez_id column - without the
                                guard, every unmapped gene would join every
                                set. Asserts all four "KEGG returned nothing"
                                causes (unresolvable organism / no pathways /
                                no Entrez conversion / no matching gene in this
                                experiment) produce their OWN note, so none of
                                them is mistaken for the generic "your DEGs
                                aren't annotated" message that would blame the
                                user's gene list for a service or coverage
                                problem. Also: progress steps_total accounts
                                for KEGG exactly, GSEA's "fewer than 2 ranked
                                genes" early-out gives KEGG its own labelled
                                block, and Reactome's output is unchanged by
                                the bridge refactor.
                                Also covers the strain-disambiguation picker
                                added on top of (1): when several name-matched
                                KEGG genomes exist and NONE taxon-verifies
                                (the real Bacillus subtilis case - 14 strain
                                genomes, none carrying species-level taxon
                                1423), resolve_kegg_organism_code now exposes
                                a structured "candidates" list (empty outside
                                this exact case) and accepts an optional
                                confirmed_organism_code. The sensitivity check
                                that matters here: a confirmed code that is
                                NOT one of THIS call's own re-derived
                                candidates is silently ignored, falling back
                                to the plain refusal - a confirmation can
                                never become an independent organism source
                                unrelated to the real taxon_name, which is
                                what would recreate the 27ce033 divergence bug
                                one layer up. Also: confirmation bypasses the
                                on-disk cache entirely (never read or
                                written), a confirmed code is a no-op outside
                                the ambiguous branch (ignored when the
                                organism already resolves cleanly), and
                                run_ora/run_gsea's new
                                kegg_confirmed_organism_code parameter reaches
                                a real result end to end - the SAME run
                                without it stays refused, proving the
                                parameter is what changed the outcome.

test_kegg_universe_parsing.py  23 checks. build_kegg_universe, the KEGG
                                PATHWAY gene-set source added to Genomics ORA
                                + Gene-set Burden Test (KEGG PATHWAY only;
                                Metabolic Connection's separate, earlier KEGG
                                work covers compound/enzyme lookups, not gene
                                sets). The universe is keyed on KEGG GENE IDs,
                                not UniProt accessions - an accession-keyed
                                first version was superseded (see
                                test_kegg_gene_map.py) because it inflated
                                counts via UniProt fan-out (measured: human
                                N=10,873 accessions vs N=9,424 real genes,
                                hsa00010 K=107 vs K=67) and could never reach
                                the ~9% of pathway genes with no UniProt
                                cross-reference at all. Guards the pathway-
                                name suffix trap found while building it:
                                names carry a " - <Organism>" suffix, but
                                several REAL pathway names ALSO contain an
                                internal " - " (e.g. "Glycosaminoglycan
                                biosynthesis - keratan sulfate") - a naive
                                rsplit(" - ", 1) mangles them, so the fix
                                strips only the single most-common suffix
                                across that organism's own rows. Also: a
                                pathway whose every gene has no link-row entry
                                is dropped entirely rather than kept empty;
                                the on-disk universe cache's "genes" key
                                round-trips list<->set correctly; and both
                                run_ora_for_database/run_burden_for_database's
                                KEGG branches produce the correct k/K/n/N and
                                result shape against a monkeypatched synthetic
                                universe (no network) - including that their
                                "Overlapping Genes"/"Genes" columns are
                                relabeled from raw KEGG gene IDs back to the
                                caller's own identifiers via kegg_id_labels,
                                falling back to the raw ID rather than
                                vanishing when no label is supplied - with
                                unknown databases still raising. Pure
                                functions and fixture strings - no server, no
                                network.

test_kegg_gene_map.py          19 checks. build_kegg_gene_map, the lookup
                                tables (UniProt accession / NCBI Entrez ID /
                                gene symbol -> KEGG gene ID) that let the
                                retriever assign each gene a KEGG gene ID -
                                the identifier build_kegg_universe's pathway
                                sets are keyed on. The one invariant a LOOKUP
                                table needs that a pathway gene-SET does not:
                                a key claimed by more than one KEGG gene is
                                dropped from the map entirely in both
                                directions (accession and Entrez bridges) and
                                for the symbol bridge - guarded by
                                _kegg_invert_unique's own sensitivity check
                                and the builder's end-to-end version.  Also:
                                _parse_kegg_gene_symbols correctly separates a
                                real "SYM1, SYM2; description" field from a
                                bare description that happens to contain a
                                ";" but carries no genuine symbols (rejects
                                space-containing and >5-token candidates); the
                                three maps are JSON-safe plain dicts (not
                                turned into sets by the shared universe cache
                                serializer); and an organism with no UniProt
                                data gets a caveat and is never cached. Pure
                                functions and a monkeypatched _kegg_safe_get -
                                no network. Uses fictitious 6+ character
                                organism codes for its fixtures and cleans up
                                any cache file it writes, specifically because
                                an earlier version of this file used the REAL
                                "cam" organism code under force_refresh=True -
                                force_refresh skips READING a stale cache, not
                                WRITING one - and silently poisoned the real
                                on-disk chickpea gene map with fake data until
                                caught by hand; this is the regression guard
                                for that incident, not a hypothetical.

test_kegg_retriever_mapping.py  20 checks. Genomics_module_retriever's KEGG
                                changes: NCBIGeneID is now kept even when the
                                UniProt bridge itself stays ambiguous (it used
                                to be thrown away, which is what made an
                                Entrez-ID-only KEGG match impossible before);
                                an ambiguous gene gets a best-effort, harvest-
                                ONLY Entrez lookup that requires exactly ONE
                                esearch hit (never gene_ids[0] of up to 5, the
                                existing not_found-only fallback's choice -
                                defensible there, but an ambiguous gene is
                                ambiguous precisely because real paralogs
                                likely exist); and resolve_kegg_ids's three-
                                bridge precedence (accession > entrez > symbol)
                                plus its accession-candidate disambiguation
                                path. The CRITICAL, headline check: a gene
                                recovered ONLY via the Entrez bridge (no
                                UniProt accession at all) must be promoted out
                                of "ambiguous" to a new status,
                                "resolved_via_kegg" - a real bug caught by
                                hand during development, where a KEGG ID got
                                filled in but the row's status silently stayed
                                "ambiguous", meaning it could never reach
                                final_df (table1_df = resolved_df[~ambiguous_
                                mask] excludes every ambiguous row
                                unconditionally) and the KEGG ID would have
                                been permanently unreachable by any analysis.
                                Also guards the inverse: two DIFFERENT genes
                                among the disambiguation candidates leaves the
                                row ambiguous, never guessing. No network -
                                _safe_uniprot_search / _ncbi_gene_search /
                                build_kegg_gene_map are all monkeypatched.

test_proteomics_retriever.py    31 checks. Proteomics_module_retriever.py -
                                the three behaviours a proteomics upload
                                needs that no transcriptomics one does:
                                protein-GROUP parsing (leading accession from
                                "P04637;P04637-2"), isoform-suffix stripping,
                                and MaxQuant contaminant/decoy filtering by
                                BOTH mechanisms at once (CON__/REV__ id
                                prefix AND "+"-flagged QC columns). Three
                                sensitivity checks, one of which reproduces
                                the SHAPE of a real bug found during
                                development - filtering a frame that never
                                carried the flag columns silently finds zero
                                flagged rows - plus a dedupe tie-break check
                                run twice with the row order REVERSED, so it
                                cannot pass by accidental first-wins. Ends
                                with a small live UniProt ID Mapping round
                                trip (real network, no local server) proving
                                a secondary accession (Q15086) resolves to
                                its current primary (P04637/TP53).

test_proteomics_ora.py          23 checks. Proteomics_module.run_ora's own
                                logic, entirely offline via a synthetic
                                Custom-database fixture (no live database
                                fetch): direction semantics, independent
                                per-direction BH families, honest empty
                                blocks, and the two statistical defects found
                                and fixed during development. BOTH are
                                sensitivity-checked by recomputing the SAME
                                fixture the pre-fix way and confirming it
                                genuinely differs: (1) zero-overlap terms
                                must never enter the BH family (unfiltered:
                                183 terms tested instead of 2); and (2) the
                                universe must be the QUANTIFIED proteins, not
                                every protein the database annotates -
                                exposed with a fixture whose terms
                                deliberately contain annotated-but-never-
                                measured accessions (pre-fix: K=25 instead of
                                the correct K=5). The original fixture was
                                blind to (2), since all of its genes were
                                already quantified. Check (1) was updated
                                2026-08-24 when Genomics_module.
                                _run_fisher_for_namespace's own zero-overlap
                                defect was fixed (b4e311f): it now asserts
                                the FIXED shared-engine behaviour (2 terms
                                tested, not 183; n_terms_no_overlap reported)
                                instead of documenting a defect Genomics
                                still had - guards against regressing it.

test_proteomics_quant.py        31 checks. Proteomics_module.
                                run_quantitative_enrichment's own logic,
                                entirely offline via a synthetic Custom-
                                database fixture with a real continuous
                                log2FC signal (not ORA's significant/non-
                                significant binary - this test has no
                                cutoff): _build_protein_scores' four metric
                                formulas, NaN-score dropping, duplicate-
                                accession mean-aggregation, and the
                                Overlapping Proteins member-column rename/
                                format shared with ORA. A real bug caught in
                                post-milestone review: n_dropped originally
                                counted dropped ROWS, not dropped PROTEINS -
                                for a duplicate accession with one valid row
                                and one NaN row (resolve_protein_table does
                                not deduplicate), it overcounted, and could
                                double-count a fully-NaN duplicate as 2
                                missing proteins instead of 1. Fixed to count
                                DISTINCT accessions with no score at all;
                                sensitivity-checked by reproducing the OLD
                                row-based formula on the same fixture and
                                confirming it genuinely gives 3, not the
                                correct 1. The LOAD-BEARING
                                check: a DOWN-shifted term must be AS
                                significant as an UP-shifted one under the
                                two-sided Wilcoxon test this milestone added
                                (Genomics_module._run_wilcoxon_for_namespace,
                                alternative param, default still "greater" so
                                the Gene-set Burden Test is unaffected) -
                                sensitivity-checked by recomputing the
                                identical fixture with the OLD one-sided
                                "greater" default and confirming the down
                                term's p-value jumps from ~1e-7 to ~1.0,
                                proving the fix is not vacuous. Also confirms
                                Genomics' own Burden Test default is
                                genuinely unchanged (omitting the new
                                `alternative` kwarg reproduces the explicit
                                "greater" call exactly).

test_proteomics_kegg.py         40 checks. KEGG in Proteomics ORA +
                                Quantitative Enrichment, entirely offline
                                (every KEGG call monkeypatched, on-disk cache
                                redirected to a temp dir). Proteomics is
                                UniProt-accession-native - unlike
                                Transcriptomics (gene-symbol-native, needed
                                an Entrez bridge), KEGG needs NO identifier
                                bridge here: "accession" is already the exact
                                key KEGG's own conv/uniprot table uses. Guards
                                build_kegg_accession_sets' accession -> gene
                                -> pathway join: the critical property is
                                that going gene -> accession reuses the SAME
                                (already-deduplicated) acc_to_gene map
                                build_kegg_gene_map produces, never KEGG's
                                raw, un-deduplicated conv/uniprot data
                                directly - that naive direction is exactly
                                the accession fan-out bug that made Genomics
                                move its OWN KEGG universe off accession
                                space (human N=10,873 accessions vs N=9,424
                                real genes). SENSITIVITY: confirms two
                                measured accessions sharing one gene (real
                                isoforms) BOTH appear - that is not fan-out,
                                it mirrors how every other Proteomics
                                accession-space database already treats
                                isoforms - while an accession dropped as
                                ambiguous by build_kegg_gene_map never
                                reappears here. Also asserts the four
                                distinct "nothing to report" notes
                                _resolve_kegg_for_proteomics (shared by both
                                run_ora and run_quantitative_enrichment) must
                                produce: organism unresolved/ambiguous, KEGG
                                unreachable, zero measured proteins mapped,
                                and mapped-but-none-on-a-pathway - each named
                                differently so none is mistaken for the
                                others or for the generic "no overlap"
                                message. Confirms the organism-resolution/
                                confirmation machinery was genuinely
                                GENERALIZED from Transcriptomics (API.
                                _kegg_confirmed_code, not a
                                Transcriptomics-only copy) and that both
                                modules' KEGG endpoints are registered.

live_proteomics_kegg.py         20 checks. LIVE network test for the same
                                feature - confirms known real human
                                accessions (GAPDH, ENO1, PKM, ALDOA -> real
                                KEGG gene IDs) land in the real
                                Glycolysis/Gluconeogenesis pathway map, and
                                CDK1/CCNB1 land in Cell cycle, through the
                                real accession -> gene -> pathway join. Also
                                re-measures the exact coverage figures the
                                module docstring cites, on the app's own real
                                3,500-protein example file: 94.6% of
                                accessions map to a unique KEGG gene, 44.1%
                                sit on a KEGG pathway map - if KEGG's own data
                                shape ever drifts enough to move these
                                numbers substantially, this is what would
                                catch it. Ends with a real run_ora +
                                run_quantitative_enrichment against live KEGG
                                membership, confirming the up-regulated
                                glycolysis genes and down-regulated cell-cycle
                                genes come out on the correct sides of both
                                engines.

live_e2e_proteomics_kegg_confirm.py
                                25 checks. Real HTTP end-to-end test (no
                                browser) for the two Proteomics KEGG
                                endpoints, mirroring
                                live_e2e_transcriptomics_kegg_confirm.py's
                                structure closely on purpose: the point is to
                                prove the GENERALIZATION (Transcriptomics-
                                only code made module-generic) didn't
                                silently narrow or break anything for the
                                second module now depending on it. Runs a
                                clean KEGG ORA + Quant pass on real human
                                accessions first, then switches to the same
                                real Bacillus subtilis ambiguity used
                                throughout this feature, confirms a real
                                candidate ("bsu") via
                                /api/proteomics/kegg/confirm-organism, and
                                proves the confirmation reaches a real run's
                                result caveat. The check that matters most:
                                re-resolves the session back to human and
                                re-runs, confirming the B. subtilis
                                confirmation does NOT leak into a run for a
                                different organism - the one layer that
                                proves API._kegg_confirmed_code's taxon
                                re-check actually guards a stale confirmation
                                at READ time, not only that
                                _job_resolve_proteins resets it at write
                                time.

test_proteomics_network.py      35 checks. Proteomics_network.py, offline
                                (fixture STRING graph, same pattern as
                                test_transcriptomics_network.py). Covers what
                                is UNIQUE to this module - the numerics
                                themselves are already covered by
                                test_network_engine_core.py and
                                test_transcriptomics_network.py, which this
                                module reuses unchanged via
                                Transcriptomics_network.run_network_analysis's
                                hook parameters (mapper/deg_selector/
                                module_enricher/entity_noun). Fixture carries
                                a PHI/PLO accession-collision pair (both
                                proteins carry "ACC_DUP"; PHI is the
                                higher-degree, correct winner - the
                                accession-space analogue of the existing
                                gene-symbol collision test) and a protein
                                (P4) reachable ONLY via the gene-symbol
                                fallback tier, never accession. Confirms
                                build_accession_string_index resolves the
                                collision deterministically and reports it;
                                map_protein_table_to_string partitions rows
                                into accession/gene_symbol tiers correctly
                                and enforces the same <10% refusal gate as
                                Transcriptomics, worded for proteins. THE
                                HEADLINE TEST: enrich_modules_by_accession
                                translates a GO term's accessions into
                                STRING-ID space BEFORE testing, so a module
                                containing PHI (which carries the colliding
                                ACC_DUP) counts k=4 (one match per PROTEIN),
                                not k=5 (one per accession) - a direct
                                sensitivity check against accession-space
                                double counting. A second test confirms a
                                RECRUITED protein (no accession in the
                                uploaded table at all) still keeps its
                                module's Reactome label, because the
                                universe is the graph, not the table -
                                mirroring Transcriptomics' own recruited-gene
                                test. Also exercises the Reactome-needs-
                                taxon_name refusal, the Complex Portal path,
                                an unsupported-database ValueError, and a
                                full run_ppi_network_analysis end-to-end
                                (24-protein synthetic graph, 12 seeds via
                                pm.select_significant_proteins, confirming
                                every seed score is strictly positive despite
                                half being down-regulated, and that the
                                returned payload carries the same _graph/
                                _hit_ids keys run_network_analysis does, so
                                the API layer needs no Proteomics-specific
                                branch to strip them). PLUS (multi-database):
                                enrich_modules_by_accession's `database` param
                                became `databases` - GO + Complex Portal
                                together returns candidates for BOTH (each
                                fixture term carries a filler, same
                                degenerate-universe fix as the Transcriptomics
                                multi-database test), Complex Portal's tighter
                                match deterministically wins while GO's own
                                candidate table (still k=4, PHI's collision
                                still resolved once) remains independently
                                usable, not discarded for having lost.

live_e2e_proteomics_network.py  31 checks. Real HTTP end-to-end test (no
                                browser) for the new /api/proteomics/network/
                                {configure,run,download} routes, mirroring
                                live_e2e_transcriptomics_network.py's
                                structure closely on purpose: the point is to
                                prove the GENERALIZATION (API._network_param_
                                bounds/_job_run_network, Transcriptomics_
                                network.run_network_analysis's hook
                                parameters) didn't silently narrow or break
                                anything for the second module now depending
                                on it. Runs the real pipeline against the
                                bundled 3,500-protein human example: real
                                UniProt ID Mapping resolve, real STRING
                                download + 500-permutation propagation, real
                                GO fetch for module labels. Confirms
                                msigdb_collections is always [] regardless of
                                what is sent (MSigDB is Transcriptomics-only);
                                an out-of-range n_permutations (200) is
                                rejected to the default (1000), never
                                honoured; an enrichment_database Proteomics
                                does not offer (MSigDB) falls back to GO
                                rather than being silently accepted; the hit
                                table carries the exact same shared column
                                set Transcriptomics' own hit table uses
                                (entity_noun changes prose only, never a
                                column name); and every CSV download
                                (hits/mapping/modules/per-module enrichment)
                                works. Mapping coverage is checked against a
                                wide floor (>50%), not a pinned historical
                                percentage - a first run measured 81.6% where
                                an earlier one-off measurement during design
                                had found ~94.7%, both comfortably above the
                                actual 10%/50% gates that matter; UniProt's
                                own live ID Mapping service is not under this
                                app's control and can vary run to run.

playwright_proteomics_network.py
                                44 checks. FIRST real-browser coverage of PPI
                                Network / Module Detection. Drives the REAL
                                flow end to end: tile -> STRING-species
                                organism selector -> upload the bundled
                                example -> resolve -> continue to PPI Network
                                / Module Detection -> configure -> run (real
                                STRING download + 500-permutation propagation
                                + GO fetch) -> results, confirming the
                                descriptor-driven generalization of
                                Transcriptomics' own Network config/results
                                screens (NETWORK_TILE_DESCRIPTORS,
                                renderNetworkConfigScreen/
                                renderNetworkResultsScreen in app.js) reaches
                                the DOM correctly rather than leaking
                                Transcriptomics wording: "Seed proteins" not
                                "Seed genes", a "Protein mapping" nav tab
                                (not "Gene mapping") whose caption reads
                                "Matched by accession...", the module-label
                                checkboxes offering Reactome AND Complex
                                Portal but NOT MSigDB, and the
                                Proteomics-specific note explaining why
                                Reactome/Complex Portal work here (unlike
                                Transcriptomics' Network tile). Also confirms
                                a real run completes, its CSV download
                                carries a "proteomics_network"-prefixed
                                filename, and there are zero console/page
                                errors across the whole flow. PLUS
                                (multi-database): checking GO AND Complex
                                Portal together does NOT uncheck either
                                (unlike the old single-select dropdown); the
                                run's own stored config echoes both databases,
                                in the order checked; every labelled module's
                                payload carries a "candidates" dict, not just
                                a single winning block; the "Module
                                enrichment" pane's own Database selector is
                                populated from a REAL run's candidates (this
                                run: GO's three namespaces plus Complex
                                Portal), switching it changes the displayed
                                annotation-source caption, and the parity fix
                                (Terms shown / Dot gradient / Bar color
                                controls, previously hardcoded with none at
                                all) landed alongside it.

test_proteomics_ptm_retriever.py
                                29 checks. Proteomics_ptm_retriever.py,
                                offline (UniProt ID Mapping monkeypatched).
                                Covers what is UNIQUE to this module - the
                                site-level input contract that
                                Proteomics_module_retriever.py cannot express
                                at all (no residue/position anywhere in that
                                module's RESOLVED_COLUMNS). normalize_site
                                handles every accepted shape (bare "S473",
                                3-letter "Ser473", phospho-prefixed "pS473",
                                lowercase, non-STY residues, split
                                residue+position columns, MaxQuant's
                                semicolon-joined multi-protein position where
                                only the first is used); detect_columns
                                accepts either a combined site column or a
                                split residue+position pair and raises
                                naming exactly what's missing otherwise. THE
                                HEADLINE TEST: isoform preservation - unlike
                                every other Proteomics tile (which strips
                                isoform suffixes via strip_isoform because
                                GO/Reactome/etc. annotate the canonical
                                entry), a phosphosite's position is
                                isoform-specific, so "P04637_S315" and
                                "P04637-2_S315" must both survive resolution
                                as DISTINCT rows rather than collapsing into
                                one; canonical_accession/canonical_site_id
                                are confirmed to be a SEPARATE pair used only
                                for external-database matching, never for row
                                identity. Also covers protein-group ambiguity
                                (a multi-member group is flagged
                                ambiguous_group=True, drop_ambiguous_groups=
                                True removes it, a single-member group is
                                not flagged), multiplicity collapse (three
                                multiplicity-class rows for the same site
                                collapse to exactly one, keeping the
                                SMALLEST-pvalue row's exact log2FC rather
                                than an average - averaging would mix two
                                distinct biological measurements),
                                localization-probability filtering (both the
                                present-and-filtering case and the
                                absent-so-honestly-False case, never silently
                                treated as satisfied), unparseable-site
                                dropping with an honest count, and reuse of
                                the existing contaminant/decoy filter.

test_proteomics_ptm.py          21 checks. Proteomics_ptm.py, offline (a
                                small synthetic iPTMnet ptm.txt fixture
                                monkeypatched in place of the real download -
                                12 real-shaped rows: KIN1 phosphorylates 6
                                sites, KIN2 phosphorylates 6 different sites,
                                plus a mouse-organism row and an
                                ACETYLATION-type row that must both be
                                excluded). Uses a RESERVED FAKE taxon id
                                ("0002"), never a real one like "9606" -
                                build_kinase_substrate_universe's parsed-JSON
                                cache (cache/ora_universe/iPTMnet_v{version}_
                                {taxon_id}.json) is SHARED with real
                                production runs and keyed only on taxon_id;
                                an earlier version of this file used the real
                                human id with force_refresh=True and silently
                                overwrote a real user's cached human kinase
                                universe (hundreds of kinases) with this
                                file's 2-kinase fixture, which then made
                                run_kinase_activity report "no kinase met the
                                minimum substrate-site count" on a
                                perfectly good real upload - bug found and
                                fixed 2026-09-21 (see also
                                playwright_proteomics_ptm.py's entry for the
                                separate chart-controls bug found in the same
                                investigation). Covers what is UNIQUE to this
                                module - the shared regulon numerics
                                (_ulm_score/_null_slopes_shared/
                                _permutation_p_value/
                                _empirical_null_calibration) are already
                                covered by test_regulon_engine.py and are
                                imported here unchanged (git diff --stat
                                Transcriptomics_regulon.py stays empty -
                                gated on explicitly), never retested.
                                build_kinase_substrate_universe confirms the
                                organism-name and ptm_type filters both work
                                (only the human PHOSPHORYLATION rows survive;
                                the mouse row and the ACETYLATION row are
                                excluded) and that a missing taxon_name or an
                                organism with zero matching rows both refuse
                                cleanly with an explanatory caveat rather
                                than raising. THE HEADLINE CLAIM, verified
                                end to end on a fixture with a KNOWN planted
                                ground truth rather than just "the code ran
                                without crashing": on a synthetic experiment
                                where KIN1's 6 substrate sites are strongly
                                up-regulated and KIN2's 6 substrate sites are
                                noise, run_kinase_activity scores KIN1 far
                                above KIN2 with a genuinely significant
                                calibrated FDR, and run_kinase_substrate_ora
                                independently finds KIN1 (only) enriched
                                among up-regulated sites via the same
                                canonical_site_id space. Also confirms the
                                min_kinase_sites gate and the "too few tested
                                sites overall" guard both return an
                                explained, empty (not None - the house
                                convention for these result tables) result
                                rather than crashing, and that per-kinase
                                regression detail is stored keyed by the
                                kinase's OWN substrate site_ids (not
                                canonical_site_ids).

live_e2e_proteomics_ptm.py      26 checks. Real HTTP end-to-end test (no
                                browser) for the new /api/proteomics/ptm/
                                {upload,resolve,configure,run,download}
                                routes, run against the bundled example file
                                (Data3/example_proteomics_phosphosites.csv -
                                300 rows built directly from REAL iPTMnet
                                data: 25 real MAPK1 substrate sites simulated
                                strongly UP, 15 real AKT1 substrate sites
                                simulated DOWN, 260 background sites as
                                noise). Runs the real pipeline: real UniProt
                                ID Mapping resolve (154/300 mapped, the rest
                                dropped mostly by the localization-
                                probability filter), real iPTMnet fetch,
                                kinase activity ULM, kinase substrate ORA.
                                Confirms isoform preservation on REAL
                                resolved data (62 isoform-suffixed sites
                                survived; spot-checked that an isoform row's
                                accession differs from, but starts with, its
                                canonical_accession). THE CENTRAL BIOLOGICAL
                                CORRECTNESS CHECK, on real data end to end:
                                MAPK1 is the top-ranked kinase in the
                                activity table (score=14.96, calibrated
                                FDR~1.5e-49, 24 substrate sites tested) and
                                AKT1 has a negative activity score (-7.07,
                                FDR~8.3e-12) - exactly matching the
                                simulated signal design; MAPK1 also tops the
                                substrate-ORA "up" table independently
                                (Fisher's FDR~1.5e-23). Also confirms
                                configure-without-resolve is refused (400),
                                all 3 CSV downloads work (activity, ora,
                                resolved), and an unknown table name is
                                rejected (400).

playwright_proteomics_ptm.py    38 checks. Real-browser coverage of the 4th
                                and final Proteomics tile, newly enabled
                                (removed from
                                PROTEOMICS_TILES_NOT_YET_IMPLEMENTED).
                                Drives the REAL flow end to end: tile
                                enabled -> UniProt-taxonomy organism search
                                (a live external lookup, unlike STRING's
                                local species file used by the other
                                Proteomics tiles, so needs a more generous
                                wait after typing the query) -> upload the
                                bundled site-level example -> resolve (real
                                UniProt ID Mapping, isoform-preserving) ->
                                continue to config -> run (real iPTMnet
                                fetch + kinase activity ULM + kinase
                                substrate ORA, 500 permutations) -> results.
                                Confirms the site-level upload section (not
                                the protein-level one) is shown for this
                                tile, MAPK1 (the kinase whose substrate sites
                                were simulated strongly UP in the bundled
                                example) is surfaced correctly in BOTH the
                                "Kinase activity" and "Kinase substrate
                                enrichment" result panes, the CSV download
                                carries a "proteomics_ptm"-prefixed filename,
                                and there are zero console/page errors across
                                the whole flow. PLUS (chart controls, bug
                                found and fixed 2026-09-21): a user reported
                                the "Kinase substrate enrichment" dot/bar
                                plots showed none of the usual chart controls
                                other tiles offer - the pane was hardcoded
                                with topN=20/fontSize=12/default colors and
                                no way to change any of them, and the kinase
                                activity bar plot likewise had no topN
                                control. Reproduced directly against the
                                real example first (isolating it from a
                                SEPARATE bug found the same day: an earlier
                                offline test run had polluted the shared
                                iPTMnet universe cache for the real human
                                taxon id with a 2-kinase test fixture - fixed
                                in test_proteomics_ptm.py by switching that
                                test to a reserved fake taxon id "0002" - see
                                that file's own header comment). Once the
                                cache was restored and the controls gap
                                confirmed as a real, separate UI bug, both
                                panes were brought up to the same feature
                                parity the PPI Network "Module enrichment"
                                pane already has (Terms shown/Rank & show
                                by/Font size/Dot gradient/Bar color). This
                                suite now confirms every one of those
                                controls is actually PRESENT on screen (not
                                just that the underlying chart-building
                                function accepts the options), that toggling
                                each one (dot gradient, font size, "Kinases
                                shown" topN, bar color) actually redraws the
                                chart without a console error, and that the
                                controls persist across a direction-key
                                switch rather than being wiped by the
                                table/chart re-render.

test_format_term_name.py        49 checks. Genomics_module.format_term_name -
                                pathway/term names rendered lowercase-first
                                ("cytoplasmic translation") across every
                                module's tables and figures, user-reported
                                2026-08-26. Root cause: GO's OWN naming
                                convention is lowercase-first (88.4% of
                                21,028 real cached GO terms); every engine
                                just passed the published name straight
                                through. A naive capitalize() is NOT safe -
                                787 real names across GO/Reactome/InterPro/
                                UniProt keywords/Complex Portal/DISEASES
                                start with a scientifically load-bearing
                                lowercase token (mRNA, cAMP, pH, eIF3,
                                hnRNP, bZIP, ...) that capitalizing would
                                corrupt into a DIFFERENT, real, wrong token
                                (cAMP -> CAMP is a different human gene; pH
                                -> PH is a different protein domain).
                                Verified against the REAL cached GO corpus,
                                not a fixture: exactly 18,228 of 21,028 real
                                terms get capitalized, exactly 352 are
                                hazard-preserved, zero false positives.
                                MSigDB's SCREAMING_SNAKE_CASE names get
                                underscores replaced with spaces, caps kept
                                (TNFA/NFKB are real gene symbols). A
                                sensitivity check proves a naive single-
                                character capitalize() - the realistic
                                accidental implementation this replaced -
                                really would corrupt 'mRNA' into 'MRNA' and
                                'pH' into 'PH'. Also integration-tests the
                                fix is actually WIRED IN, not just correct
                                in isolation: calls Genomics_module.
                                _run_fisher_for_namespace (ORA) and
                                _run_wilcoxon_for_namespace (Burden Test /
                                Proteomics Quantitative Enrichment, which
                                reuses it) and Transcriptomics_module.
                                _run_ora_for_sets directly and confirms both
                                a normal name and a hazard name come out
                                correctly in the real result table.

test_metabolomics_concentration.py
                                58 checks. Metabolomics_module's concentration-
                                table pipeline (parse_concentration_table,
                                compute_compound_statistics,
                                build_scores_from_stats), offline. Exists
                                because the first Metabolomics input contract
                                was WRONG: it asked for a user-supplied
                                log2FC/pvalue/padj per compound (a
                                differential-expression convention metabolomics
                                does not use). Covers both table orientations
                                (samples-in-rows and the transposed
                                compounds-in-rows lipidomics layout), group-
                                column detection and its honest failures (one
                                group, no label column), explicit overrides
                                being honoured or rejected - never silently
                                replaced by the auto-detected column - a
                                planted 4x shift recovered as log2fc > 1.5 with
                                p < 0.001 while an unchanged compound is not
                                significant, min_per_group dropping compounds
                                (counted, never imputed), and the score handed
                                to the Wilcoxon engine matching the computed
                                statistic.
test_metabolomics_crosswalk.py  29 checks. The identifier-bridging layer,
                                offline (bulk fetches monkeypatched to a small
                                fixture): source merge and graceful loss of an
                                unreachable source, the match-tier order,
                                tier counts partitioning matched rows, the
                                InChIKey-skeleton fallback being counted
                                separately and never overriding an exact
                                match, and lipid shorthand being class-yes /
                                pathway-no. PLUS the universe_scope split:
                                experiment scope narrows to the user's own
                                compounds, database scope covers the whole
                                crosswalk, and a lipid SHORTHAND compound is
                                unioned into its own class under database
                                scope (an earlier version selected these by a
                                "SHORTHAND:" key prefix, silently dropping a
                                shorthand name that RefMet also lists - caught
                                by this file).
test_metabolomics_ora.py        62 checks. run_ora (compound-LIST input,
                                database universe) and
                                run_quantitative_enrichment (measured
                                universe), offline. Only the KEGG download and
                                the universe cache are faked, so the REAL
                                builder and overview-map exclusion run (an
                                earlier version replaced the builder with a
                                fake that re-implemented the exclusion, which
                                tested the fake). Covers: KEGG map01000-01099
                                never entering a result table and the real TCA
                                pathway winning; the overlap column showing
                                the user's own names, semicolon-joined, with a
                                comma-containing name surviving as ONE member;
                                no up/down split; database universe larger
                                than measured; a background list narrowing it;
                                Chemical Class from the crosswalk taxonomy.
                                KEGG is ORGANISM-SPECIFIC (KEGG serves no
                                organism-specific compound links, and a map's
                                compounds do not vary by organism, so the
                                organism selects which MAPS are tested): terms
                                carry the organism's ids and names, a map the
                                organism lacks is excluded and shrinks the
                                universe (asserted on the universe, not just
                                the rows), a second organism gets a DIFFERENT
                                term set from the SAME shared link data, the
                                overview-map exclusion still holds AFTER
                                re-keying to hsa..., no genome / ambiguous
                                organism -> EMPTY with a note and no fallback,
                                a confirmed code reaches the resolver, and the
                                organism caveat survives beside a note. An
                                EXPLICIT KEGG organism (the KEGG search box) is
                                used as given with the organism derivation not
                                consulted at all (a stub that raises proves
                                it), wins over an ambiguous derived organism
                                through run_ora AND run_quantitative_
                                enrichment, and a code KEGG has no maps for is
                                empty with its own note; the failure notes are
                                re-worded for Metabolomics (the shared
                                resolver's say "try GO or MSigDB"). The
                                resolver itself is stubbed (it reaches KEGG
                                through Transcriptomics_module's own alias,
                                which faking this module's _kegg_safe_get does
                                not intercept). Two deliberate mutants - the
                                organism filter removed, and _empty_block
                                dropping the caveat - were confirmed to turn
                                these checks red. Organism ids are read from the
                                5-digit TAIL (4-character codes, 'path:' prefix).
test_metabolomics_retriever.py  14 of 15 checks pass. Metabolomics_retriever.py
                                (a deliberate independent fork of
                                test_compound_resolver.py so Metabolic
                                Connection can never be disturbed): its own
                                cache directory (cache/metabolomics), tolerance
                                of empty/corrupt cache files, the output
                                column contract. KNOWN FAILURE, not caused by
                                the fork: "HMDB IDs can be enriched from
                                PubChem synonym data" fails identically in the
                                ORIGINAL test_compound_resolver.py - the
                                inherited inline test is not isolated from
                                whatever real cache state is on disk (the CID
                                it mocks resolves from cache instead).
test_metabolomics_session.py    16 checks. Direct in-process tests of API.py's
                                job functions (no live server, no network -
                                same pattern as test_organism_consistency.py),
                                covering two session-state bugs found after
                                the module shipped. FIX 1: the LIST tiles
                                (ORA/Chemical Class/Lipid Class) never
                                deduplicated an exact-repeated compound Name -
                                remove_duplicate_input_names existed but was
                                only ever called from the CLI resolver, never
                                from the API path. Now wired in for the list
                                path only (dedupe_names=True by default);
                                confirmed the CONCENTRATION path still never
                                dedupes, since a repeated column Name there is
                                a genuine second feature/adduct, not a
                                duplicate. FIX 2: the list and concentration
                                paths share ONE organism (taxon_id/taxon_name),
                                but each resolve job only cleared its OWN
                                downstream results - resolving one path for a
                                new organism left the other path's cached,
                                organism-dependent results silently stale.
                                Now resolving either path for a genuinely NEW
                                organism also clears the other path's
                                organism-dependent results (ORA <-> Quant),
                                while re-resolving for the SAME organism
                                leaves them alone, Chemical/Lipid Class are
                                left untouched either way (they never use the
                                organism), and quant_stats_df/meta (per-
                                compound statistics) are never cleared by an
                                organism change. Both fixes were mutation-
                                tested: reverting either one in API.py turns
                                the corresponding checks here red.
test_metabolomics_topology.py  43 checks. Pathway Topology Analysis
                                (Metabolomics_module.run_pathway_topology),
                                offline - only KEGG's network edges are
                                faked, so the REAL enrichment path
                                (build_kegg_compound_universe, _eligible_
                                keys, _drop_zero_overlap, _run_fisher_for_
                                namespace / _run_wilcoxon_for_namespace) runs
                                unchanged, reused from test_metabolomics_
                                ora.py's own KEGG coverage. HEADLINE
                                GUARANTEE: run_pathway_topology calls the
                                exact same building blocks run_ora / run_
                                quantitative_enrichment use for their own
                                KEGG entry, so p-values/FDRs are checked
                                BYTE-IDENTICAL against a live run_ora /
                                run_quantitative_enrichment call on the same
                                fixture, in both input modes (list -> Fisher,
                                concentration table -> Wilcoxon). Also:
                                parse_kegg_kgml (reversible reactions add
                                both edge directions, irreversible only one,
                                a multi-id KGML entry becomes one node per
                                id, glycan gl: ids are kept as graph nodes,
                                a real map with zero <reaction> elements
                                returns an empty edge list rather than an
                                error); _brandes_betweenness (this module's
                                own implementation, no networkx) checked
                                against a BRUTE-FORCE shortest-path
                                enumeration on 25 random directed graphs;
                                pathway_impact's well-defined flag
                                (SENSITIVITY: a zero-reaction map's impact is
                                0.0 but explicitly flagged, never silently
                                indistinguishable from a genuine zero);
                                fetch_pathway_graph's cache (a 'nodata' 404
                                is cached, an 'error' never is) and its
                                organism-specificity (SENSITIVITY: the SAME
                                map number for two different organisms
                                yields two DIFFERENT reaction graphs, proving
                                the cache key includes the organism); and a
                                no-KEGG-genome organism returning empty with
                                a note rather than a fallback. Also:
                                build_kegg_compound_name_index / fetch_
                                pathway_graph label each node with its REAL
                                compound name (from KEGG's own list/compound
                                bulk file, cached separately from the KGML
                                fetches) rather than the bare C-number, with a
                                bare-id fallback for a node absent from the
                                index; the fake KEGG responder's /list/
                                compound handling is excluded from the KGML
                                fetch-count assertions above (a separate,
                                independently-cached artifact).

test_multiomics_joint.py        96 checks. Multi-omics Joint Pathway Enrichment engine
                                (Multiomics_module), offline, fixtures only.
                                HEADLINE GUARANTEE: the per-omics p-value inside a
                                joint result equals the one that omics' OWN ORA
                                engine produces for the same input, checked
                                against the exact function each module uses
                                internally (Genomics' _run_fisher_for_namespace,
                                Transcriptomics' _run_ora_for_sets, Metabolomics'/
                                Proteomics' _drop_zero_overlap + _run_fisher_for_
                                namespace). The design rules, each mutation-tested
                                (a mutant breaks the named check, then the file is
                                restored byte-identical): (1) a zero-hit pathway
                                keeps p = 1.0 for the combination, whereas the
                                single-omics engine drops it; (2) a pathway nobody
                                hit never enters the combined family; (3) a
                                pathway tested in only some omics goes to a partial
                                table and is never padded with p = 1.0; (4) an
                                omics with no annotated hit is left out and fewer
                                than two usable omics gives an empty result with a
                                note; (5) Benjamini-Hochberg is applied ONCE, over
                                the combined p-values only (a mutant applying it
                                twice fails); (6) Stouffer is not offered - the
                                test shows WHY (one omics at p = 1e-10 plus one at
                                p = 1.0 gives p = 0.90 under Stouffer, 2.4e-9 under
                                Fisher). Fisher equals a hand-computed -2*sum(ln p)
                                against chi-squared on 2k df; Tippett equals
                                1-(1-min p)^k. BROWN'S METHOD (the default; Brown
                                1975, empirical version Poole 2016), each property
                                checked and mutation-tested (7 engine mutants, all
                                caught): rho = 0 equals Fisher exactly; an omics
                                combined with an exact copy of itself returns its
                                own p (Fisher does not); Brown >= Fisher wherever
                                Fisher's p <= 0.28 (an honesty test: above ~0.28
                                the tails cross, so the claim is NOT made there);
                                negative rho clamped to 0; rho = Pearson of
                                -2 ln p over ALL jointly tested pathways, checked
                                against a pure-python Pearson (numpy.corrcoef
                                crashes the interpreter in this environment); fewer
                                than 20 pathways falls back to Fisher with a note.
                                REGRESSION: Transcriptomics' symbol
                                sets carry 'genes' as a sorted LIST, not a set - the
                                first draft crashed on it. The per-omics input
                                builders are checked to call that omics' own
                                universe builders and to pass the ONE shared KEGG
                                organism code through.

test_multiomics_mapping.py       124 checks. Cross-omics Pathway Mapping engine
                                (Multiomics_mapping) and the structure fetchers in
                                Multiomics_module_retriever, offline - small embedded
                                KGML / Reactome diagram+graph fixtures and mocked HTTP.
                                The rules, each mutation-tested (38 mutants across this
                                file, the retriever and API.py: all caught, files
                                restored byte-identical): the KGML parser keeps GENE
                                entries and the reaction -> enzyme link (Metabolomics'
                                parser drops both), a glycan is never a compound, a
                                multi-gene box keeps every gene; the Reactome parser
                                opens a complex / set to its UniProt / ChEBI leaves,
                                guards cycles, strips an isoform suffix, keeps a
                                cofactor inside a complex OUT of the metabolite members,
                                treats a set of only small molecules as a metabolite,
                                and restricts a sub-pathway to its own reactions;
                                three states per omics (hit / measured / absent, never
                                two) and a box is hit if ANY member is; every match
                                records its identifier bridge (Entrez / UniProt / gene
                                symbol -> KEGG gene, KEGG gene -> UniProt, symbol only
                                when no UniProt is known, ChEBI, and the InChIKey
                                skeleton for another protonation state - which never
                                overrides an exact match and never matches a different
                                skeleton); a reaction link needs enzyme-hit +
                                metabolite-hit or enzyme-hit in two omics, and a
                                'concordant' / 'discordant' call needs both expression
                                omics; the picker ranks by omics count FIRST and drops
                                zero-hit pathways; the fetchers never cache a transport
                                error, do cache a definite 'no map', and draw a pathway
                                with no diagram from its NEAREST ancestor.
                                SECOND ROUND (the maps the user opened and found
                                wrong): KEGG's global overview maps draw a gene as a
                                LINE with no x / y - the parser keeps the polyline
                                (corner points de-duplicated) and the middle vertex, a
                                line ortholog is decoration, never a node; a group (a
                                drawn complex) holds the union of its components' genes
                                and is hit if any is; KGML relations are read with every
                                subtype (a relation to a missing entry and a maplink
                                are skipped; 'compound' is not a type); interaction
                                links need BOTH ends hit; a drawn complex is not a box
                                and a gene drawn in several boxes is ONE member hit;
                                Reactome overview pathways (only ProcessNodes) are
                                detected, their sub-pathways resolved one by one, an
                                unreachable lookup is an error and NOT cached, a
                                sub-pathway reached twice is listed once. 22 more
                                engine / retriever mutants, all caught.
                                THIRD ROUND (readable big maps): a line gene keeps
                                KEGG's colour for it (upper-cased: '#b3b3e6' and
                                '#B3B3E6' are one colour) and a pathway label the map
                                number it links to; KEGG's classification (br08901) is
                                read at its B level for the METABOLISM branch only
                                (disease pathways share metabolic enzymes and once named
                                a colour 'Neurodegenerative disease'); a colour is named
                                by the class most of its genes' pathways belong to, the
                                overview maps do not vote, a tie goes to the
                                alphabetically first class, an unclassified colour is
                                'Other'; the legend lists each class once with its lines
                                and hits and its most frequent colour as the default. 9
                                more mutants, all caught.
test_multiomics_mapping_session.py  69 checks. The mapping routes and jobs in API.py,
                                offline (no server, no network). State and clearing (a
                                re-resolve, a re-upload and a re-configure drop the
                                mapping); validation (MetaCyc / no database, a
                                malformed pathway id - '../../etc/passwd', a Reactome
                                id for KEGG - is refused before it can reach a URL, and
                                no job starts); the picker and the run use the SAME
                                hit rule as Joint (shared helpers), resolve the KEGG
                                organism ONCE, skip Genomics' KEGG-id pass for a
                                Reactome map, relabel KEGG compounds, report a missing
                                organism / a missing pathway / a transport error
                                explicitly, are strict-JSON safe, and leave the four
                                single-omics session sub-dicts untouched. Also: the
                                payload carries only the omics that touched a box
                                ('absent' is left out - a global map has ~5,000 boxes),
                                a line gene's polyline and a group's components, the
                                relations and the interaction links (and a CSV for
                                them; an unknown table name never reaches a path); an
                                overview result lists its sub-pathways with the
                                picker's hit counts, building the counts itself when no
                                picker exists yet. 7 more API mutants, all caught.
                                THIRD ROUND: a line map's result lists its pathway
                                classes (KEGG colour, lines, hits), every line gene and
                                every metabolic pathway label carries its class, a map of
                                boxes and a Reactome diagram have none. 3 more mutants.
test_multiomics_retriever.py    22 checks. Multiomics_module_retriever (the thin glue file),
                                offline. It must add NO resolution logic of its
                                own: every upload / resolve call is routed, with
                                spies on the real modules, to that omics' OWN
                                loader / resolver with the arguments its own tile
                                passes (Genomics gene symbols get the organism,
                                Ensembl IDs do not; Transcriptomics auto-detects
                                its identifier type). DE selection pools up AND
                                down and equals Transcriptomics' select_degs /
                                Proteomics' select_significant_proteins for
                                direction 'all'; the reference set is everything
                                TESTED, not a database. Ambiguous / unresolved
                                Genomics rows are dropped, never guessed.
                                Mutation-tested: a direction-split selection,
                                keeping ambiguous genes, a forced Transcriptomics
                                identifier type, and a significant-only reference
                                set each turn a named check red.

test_multiomics_session.py      47 checks. The Multi-omics session state and jobs in API.py,
                                no server: route functions and jobs are called
                                in-process with stubs. ISOLATION (the user's
                                explicit concern): a full select -> resolve ->
                                configure -> run pass leaves session['ora'] /
                                ['transcriptomics'] / ['proteomics'] /
                                ['metabolomics'] untouched. ONE KEGG ORGANISM per
                                run: resolve_kegg_organism_code is called exactly
                                once for a four-omics, two-database run and the
                                SAME code reaches every omics' KEGG builder and
                                Genomics' KEGG-ID pass; with no usable organism no
                                KEGG builder is called at all (no fallback) while
                                Reactome still runs; a stale confirmation made for
                                another organism is ignored. A failed resolve
                                commits nothing (organism, earlier omics' tables
                                and results all stay as they were). Route
                                validation: <2 omics, out-of-order steps, Stouffer,
                                bad cutoffs / term sizes. Mutation-tested: a
                                per-omics organism resolution, committing the
                                organism before the resolvers run, KEGG builders
                                still called without an organism, and the run job
                                writing into the Metabolomics session each fail a
                                named check.


---------------------------------------------------------------------------
NEEDS A LIVE SERVER on http://127.0.0.1:8009
---------------------------------------------------------------------------
Start it first (from the repo root), and note there is deliberately no
--reload, so RESTART IT after any Python edit:

    python -m uvicorn API:app --port 8009

smoke_regulon_fetchers.py       21 checks. Real TRRUST / OmniPath / ChEA3
                                fetches. Uses a throwaway cache dir so it
                                cannot poison the real one. Also an
                                UPSTREAM-COVERAGE WATCH: probes realistic
                                next-candidate organism names (rat,
                                zebrafish, fly, ...) against TRRUST's real
                                URL pattern and confirms they still 404 -
                                TRRUST_URL_BY_TAXON hardcodes exactly
                                human+mouse because that is genuinely all
                                TRRUST has ever published; the moment one of
                                these stops 404ing is the signal to widen the
                                gate, not something that should go unnoticed.
                                A sensitivity check runs the identical
                                detection logic against "human" (which
                                genuinely exists, HTTP 200) and confirms it
                                is correctly flagged non-404, using real live
                                data rather than a fabricated response.
live_e2e_regulon.py             15 checks. Full job lifecycle through the API.
live_e2e_regulon_analyses.py    12 checks. Enrichment-only, activity-only and
                                both-analyses payload shapes.
live_e2e_example_ensembl.py     12 checks. The Ensembl example file end to end.
regression_ora_gsea.py          22 checks. The pre-existing GO-only ORA/GSEA
                                paths still work, PLUS real HTTP-level
                                coverage of MSigDB through the actual API:
                                configure validation ("MSigDB" without
                                msigdb_collections is rejected), a real job
                                run with a real MSigDB download, the new
                                "MSigDB: <collection> (direction)" result-key
                                convention, and n_independent_families in the
                                summary. A different layer than
                                test_msigdb_fetcher.py (direct function
                                calls) or playwright_regulon.py (browser).

live_e2e_transcriptomics_network.py  25 checks. Full real HTTP round trip
                                for Transcriptomics Network / Module
                                Detection - upload -> resolve -> configure ->
                                run -> download x3 -> guards - against the
                                bundled human example, at 500 permutations
                                (200 is now confirmed REJECTED to the new
                                default of 1000, not silently honored - a
                                dedicated check for this). Real mygene.info/
                                STRING/QuickGO calls. Under the corrected,
                                production-scale-calibrated statistics (see
                                calibration_network_production_scale.py),
                                this specific example legitimately produces 0
                                significant proteins at 500+ permutations -
                                assertions on hit/recruited/module counts and
                                caveat presence are REPORT-ONLY or gated on
                                the real observed count rather than assuming
                                a positive result, so the suite stays honest
                                either way. Confirms mapping coverage,
                                transcriptomics-shaped hit columns (no
                                "UniProt ID"), a real subnetwork edge list,
                                direction-composition columns, full
                                reproducibility parameters, and each of the
                                (up to) 3 CSV downloads working; guards an
                                unknown table name, a missing resolved table,
                                MSigDB labelling with no collection selected,
                                and (new) that a "which=modules" download
                                correctly returns 400 when a run produced no
                                modules. live_e2e_transcriptomics_network_
                                mouse.py below is the reference for what a
                                genuine positive result looks like under this
                                same, fully-corrected engine.

live_exercise_transcriptomics_network_branches.py  16 checks. Exercises
                                branches of the Network module that existed
                                but had never actually RUN even once: a real
                                abs_metric-weighted run (a real winsorization
                                quantile, weight_ratio<=20), real
                                direction=up/down runs (each a genuine SUBSET
                                of the direction=all seed count, proving
                                direction restricts selection rather than
                                relabelling it), a real MSigDB-labelled
                                module run, the enrichment CSV download
                                (called here for the first time) plus its
                                unknown-module-id guard, and the Ensembl/
                                Entrez identifier fallback tiers - TP53 and
                                MYC's REAL STRING protein ids/aliases, read
                                directly from the already-cached human
                                aliases file, paired with gene_symbol values
                                STRING can never recognise, confirming
                                map_de_table_to_string resolves them ONLY
                                through the alias tiers and lands on the
                                correct protein. Uses n_permutations=500 (was
                                100, itself always an invalid value silently
                                rejected to whichever default was active -
                                corrected here to the real evidence-supported
                                floor, restoring the "keep it fast" intent
                                without asserting an invalid setting).

live_e2e_transcriptomics_network_mouse.py  12 checks. A full real run on a
                                SECOND organism (mouse, taxon 10090), kept in
                                its own file so a mouse STRING download
                                failure is diagnosable on its own rather than
                                aborting the human-only suite - the first
                                real, non-human exercise of this module end
                                to end. Data3/example_transcriptomics_
                                mouse_test.csv (30 well-known mouse genes,
                                real symbols, synthetic log2FC/p-values).
                                Confirms taxon 10090 actually drives the run
                                (not a leftover human session), a real
                                ~16,000-protein mouse network gets built,
                                29/30 genes map, and hits/modules are
                                produced (13 hits, 12 recruited at time of
                                writing, under the same corrected, skew-aware
                                statistics as the human example above) - the
                                canonical proof that this module's positive
                                capability still works end to end, now that
                                the bundled human example alone legitimately
                                shows 0 hits at production-calibrated
                                permutation counts.

live_msigdb_urls.py             49 checks. LIVE network test: every pinned
                                MSigDB URL (14 available collection/organism
                                pairs) actually resolves, and its gene-set
                                count matches MSigDB's own published figures
                                (or, where none exists, this session's own
                                directly-measured baseline) within a stated
                                tolerance. A sensitivity check deliberately
                                breaks the version pin and confirms the fetch
                                fails loudly - this is what would have caught
                                the OLD Hallmark-via-Enrichr fetcher's silent
                                6-year staleness. Also an UPSTREAM-COVERAGE
                                WATCH: fetches MSigDB's own release index for
                                the current pinned version and confirms it
                                still lists exactly {Hs, Mm} and nothing else
                                - MSIGDB_TAXON_SUFFIX hardcodes only those two
                                because that is genuinely all MSigDB has ever
                                published (verified live across every release
                                2022.1 through the current pin); if a third
                                organism ever appears, this goes red rather
                                than silently leaving the app narrower than
                                upstream. A synthetic-HTML sensitivity check
                                proves the parsing/comparison logic itself
                                would catch a third organism, since this
                                suite cannot control what the real remote
                                server publishes. Does NOT need the app's own
                                uvicorn server - only real internet access to
                                MSigDB's release server, same as
                                smoke_regulon_fetchers.py below.

live_transcriptomics_kegg.py   31 checks. LIVE network test for the same
                                feature - what the offline file cannot prove:
                                that its fixtures match what KEGG actually
                                serves. Resolves 9 real organisms to their
                                KEGG codes and taxon-verifies each (human,
                                mouse, rat, zebrafish, fly, worm, Arabidopsis,
                                chickpea, and yeast via the genus+species
                                fallback), confirms the real Bacillus subtilis
                                ambiguity (14 strain genomes, none carrying
                                taxon 1423) refuses rather than guesses, and
                                measures the Entrez bridge on real data: 372/372
                                human pathways survive it, chickpea has 163
                                pathways over 5,757 genes (the non-model case
                                KEGG is offered for), and every member ID comes
                                out as a bare Entrez ID. Ends with a real
                                run_ora + run_gsea against live KEGG membership.
                                Also confirms the strain-disambiguation picker
                                against the real B. subtilis candidate list:
                                its "candidates" carry all 14 real strains
                                (not just the note text), confirming the real
                                "bsu" fetches a genuine KEGG pathway set
                                (128 pathways) exactly like any automatically-
                                resolved organism, and confirming a code KEGG
                                never actually name-matched ("hsa") is ignored
                                against the LIVE candidate list too, not just
                                the offline fixture.

live_e2e_transcriptomics_kegg_confirm.py
                                18 checks. Real HTTP end-to-end test (no
                                browser) for the two new endpoints the picker
                                calls - /api/transcriptomics/kegg/resolve-
                                organism and /confirm-organism - mirroring
                                regression_ora_gsea.py's approach. Uses the
                                same real Bacillus subtilis ambiguity as
                                live_transcriptomics_kegg.py. The check that
                                matters most here: confirms a KEGG organism
                                for B. subtilis, runs ORA over HTTP and sees
                                the confirmation's caveat in the result: THEN
                                re-resolves the session to a different
                                organism (human) and re-runs - proving the
                                confirmation does NOT leak into a run for a
                                different species. This is the one layer
                                neither the offline nor the direct-call live
                                test can see: whether API.py's
                                _transcriptomics_kegg_confirmed_code actually
                                guards a stale confirmation at read time, not
                                only whether _job_resolve_de resets it at
                                write time. Also: resolve-organism is
                                confirmed read-only (identical result on two
                                consecutive calls), and confirm-organism
                                rejects both a non-candidate code and a
                                missing one.

live_kegg_universe.py           39 checks. LIVE network test for KEGG
                                PATHWAY in Genomics ORA + Gene-set Burden
                                Test: build_kegg_universe, build_kegg_gene_map
                                and Genomics_module_retriever.resolve_kegg_ids
                                together, against real KEGG/UniProt/NCBI data.
                                Confirms the universe's re-keying onto real
                                KEGG genes (human N=9,424, hsa00010 K=67 -
                                DELIBERATELY down from an earlier accession-
                                keyed version's inflated 10,873/107, not a
                                regression - see test_kegg_gene_map.py for
                                why); chickpea (non-model organism) coverage;
                                a real total-failure organism (list/pathway
                                200-empty AND link/pathway HTTP 400); the
                                measured symbol-coverage asymmetry that ruled
                                out a GO-style symbol fallback as the primary
                                fix (human ~94% of genes carry a KEGG symbol,
                                chickpea ~1.9%); a REAL organism with pathway
                                data but no NCBI Entrez cross-reference
                                (Synechocystis/"syn" - has_entrez=False, still
                                builds and caches); and the identity
                                assumption this whole feature leans on (a KEGG
                                gene ID's numeric part equals its Entrez ID)
                                checked against 200 real chickpea genes, not
                                assumed. The HEADLINE check: the exact real-
                                world case that motivated this feature - real
                                chickpea gene symbols ARF2/ARF6 come back
                                "ambiguous" under UniProt's own reference-
                                proteome tie-break (no reference proteome
                                exists for most non-model organisms) - run
                                through resolve_gene_table then
                                resolve_kegg_ids and confirmed EVERY one of 7
                                real genes ends up with a KEGG gene ID, zero
                                left ambiguous. A separate KNOWN-ANSWER check
                                sends 10 real human gene SYMBOLS through the
                                REAL retriever pipeline (not hand-typed
                                accessions) into the app's REAL Fisher engine
                                and confirms hsa00010 "Glycolysis /
                                Gluconeogenesis" ranks 1 (FDR ~5e-21, all 10
                                matched, zero dropped, K=67 - the real KEGG
                                gene count), with "Overlapping Genes" showing
                                real gene symbols via kegg_id_labels, never a
                                raw KEGG gene ID, and the "Metabolic pathways"
                                overview map excluded by max_term_size rather
                                than a build-time exclusion. Cross-checked
                                against an INDEPENDENT database in a DIFFERENT
                                ID space - the same genes, resolved via
                                UniProt, also top Reactome's own Glycolysis
                                pathway - and against the Burden Test engine
                                (with its own "Genes" column relabeling
                                confirmed) on the same known input. Does NOT
                                need the app's own uvicorn server - only real
                                internet access to KEGG's REST API
                                (rest.kegg.jp) and UniProt/NCBI (via the real
                                retriever), same pattern as live_msigdb_urls.py
                                above.

live_e2e_metabolomics.py        97 checks. Real HTTP end-to-end test (no
                                browser) of all five Metabolomics tiles:
                                upload -> resolve (real PubChem/KEGG name
                                resolution + the local LIPID MAPS/RefMet
                                crosswalk) -> configure -> run -> download.
                                Asserts the BIOLOGY the bundled examples plant
                                rather than shapes: ORA's KEGG top term is a
                                real pathway (not an overview map) against the
                                database universe, Chemical Class recovers
                                "TCA acids", Lipid Class reaches
                                Glycerophosphocholines from shorthand PC
                                names, and Quantitative Enrichment finds the
                                planted TCA/glycolysis increase UP and the lipid
                                decrease DOWN (two-sided). Also the input-
                                contract errors that motivated the rework: a
                                list with no Name column, a single-group
                                table, an unknown orientation and an unknown
                                group column each rejected with a message.
                                KEGG organism-specificity against REAL KEGG:
                                human terms are hsa-prefixed with a caveat
                                naming the organism, yeast gets a strictly
                                smaller universe from the same shared link
                                file, an organism KEGG has no genome for is
                                EMPTY with a note (never a fallback), and B.
                                subtilis (14 KEGG strains) is empty until a
                                KEGG organism is picked explicitly:
                                /api/metabolomics/kegg/check-organism returns
                                the organism's name and map count and rejects an
                                unknown or malformed code, configure
                                re-validates it, an explicit pick (bsu) makes
                                KEGG test bsu maps, and beside a human
                                organism an explicit yeast pick still tests
                                yeast maps with its own smaller universe.
                                Also cross-tile staleness: resolving Quant for
                                a different organism (mouse) makes ORA's
                                already-computed KEGG results undownloadable
                                (400) until ORA is re-run for the current
                                organism, and a repeated compound name in an
                                uploaded list is kept at upload time but
                                removed exactly once at resolve time, with
                                n_input_rows still reporting the original
                                (pre-dedup) count. Pathway Topology Analysis
                                (KEGG only, both input modes): checked
                                directly against a fresh ORA / Quantitative
                                Enrichment KEGG run on the SAME input -
                                p-values and FDRs must be BYTE-IDENTICAL
                                (run_pathway_topology calls the exact same
                                building blocks internally), every row
                                carries Impact in [0, 1], a graph is fetched
                                for every tested pathway and only those, and
                                a no-KEGG-genome organism is empty with a
                                note. Needs network.

live_e2e_multiomics_mapping.py  93 checks. Real HTTP end-to-end test of Cross-omics Pathway
                                Mapping (port 8009, real network, fresh session per
                                run): four omics uploaded and resolved, the KEGG and
                                Reactome pickers, KEGG hsa00020, Reactome R-HSA-71403
                                and a sub-pathway with no diagram of its own
                                (R-HSA-70171, drawn inside Glucose metabolism). It
                                recomputes INDEPENDENTLY from the raw files: the
                                Genomics hits on the KEGG map equal the genes found by
                                re-reading KEGG's KGML and gene list (18 of 29), and
                                every Transcriptomics member's log2FC and hit/not-hit
                                equals the uploaded CSV row at the cutoffs. Also: a
                                stricter cutoff never adds a hit and turns dropped
                                genes into 'measured', not 'absent'; citrate synthase
                                is linked to citrate; no protein / complex box is ever
                                a Metabolomics hit; MetaCyc and malformed ids are
                                refused; a re-upload clears the mapping. Two of its
                                first assertions were WRONG, not the code - the bundled
                                expression tables plant no TCA signal, so TCA is hit by
                                Genomics + Metabolomics (2 omics, rank 69 of 340 as the
                                ranking puts omics count first) - and were corrected to
                                what the data supports rather than loosened.
                                SECOND ROUND, on the maps that looked wrong: hsa01100
                                (all 1,342 genes carry a polyline starting where KEGG's
                                own coords start - recomputed from KEGG's raw file - and
                                none sits at (0, 0); the Genomics hits equal a re-read of
                                the file; the payload is ~1.8 MB), hsa04210 (all 135
                                relations and the 4 groups match KEGG's file; the 14
                                interaction links are exactly the relations with both
                                ends hit, recounted from the payload) and Reactome's
                                Disease / Metabolism / Immune System (reported as
                                overviews with their sub-pathways, whose counts equal
                                the picker's; Reactome nests overviews, so the test
                                follows them down - Infectious disease > Viral
                                Infection Pathways > Dengue ... - to a real diagram).
                                THIRD ROUND: hsa01100 has 9 pathway classes, every one a
                                class of KEGG's metabolism classification read
                                INDEPENDENTLY from br08901; all 1,342 line genes carry
                                one and each of KEGG's 14 line colours maps to exactly
                                one class; citrate synthase (a TCA-cycle enzyme) comes
                                out as Carbohydrate metabolism; all 169 pathway labels
                                carry the class br08901 gives their map.
live_e2e_multiomics.py          126 checks. Real HTTP end-to-end test of Multi-omics Joint
                                Pathway Enrichment (no browser): choose omics ->
                                upload one file each -> resolve with each omics'
                                OWN resolver (UniProt / mygene.info / UniProt ID
                                mapping / PubChem + crosswalk) -> KEGG organism ->
                                configure -> run -> download. Asserts the BIOLOGY
                                the examples plant (TCA-cycle genes + TCA-cycle
                                compounds): KEGG's TCA cycle and glycolysis are
                                jointly significant with a combined p-value
                                stronger than either omics alone (TCA: 7e-50 vs
                                5.6e-41 for the better single omics), and Reactome
                                puts a TCA / glycolysis pathway in its top results.
                                HEADLINE: in a SECOND session the real Metabolomics
                                ORA tile (KEGG + Reactome) and the real Genomics
                                ORA tile (KEGG + Reactome) are run on the same
                                input, and every per-omics p-value in the joint
                                result equals theirs (46 pathways compared). Also:
                                Tippett is a real second method (1-(1-min p)^2),
                                Stouffer is refused, a four-omics run with the
                                Transcriptomics and Proteomics DE tables (cutoffs
                                applied by Omiscera, up and down pooled; 567 of
                                2969 genes significant) combines four p-values per
                                row, a plain gene list is REJECTED for
                                Transcriptomics, the partial-coverage table never
                                carries a combined p-value, one KEGG organism is
                                reported for the whole run, and the downloads are
                                real CSVs. Every p-value comparison uses a RELATIVE
                                tolerance: an earlier absolute one (1e-9) could not
                                tell 1e-40 from 1e-30, so the tile-equality checks
                                were vacuous for strong pathways. Each run uses a
                                fresh session id (the server keeps sessions, and
                                re-choosing omics keeps their uploads). Needs network.

---------------------------------------------------------------------------
NEEDS A LIVE SERVER **AND** PLAYWRIGHT (real Chromium, real clicks)
---------------------------------------------------------------------------
    pip install playwright && playwright install chromium

playwright_regulon.py           90 checks. Drives the actual UI: retriever ->
                                config -> results, asserting each figure's
                                traces and encodings, that the results nav
                                sits inside the standard '.ora-db-card' box
                                (blue top accent) matching GSEA/ORA - a real,
                                user-reported bug: this screen used to append
                                its nav straight to the page with no wrapping
                                card, visually inconsistent with every other
                                results screen - that the Regulon Detail
                                pane's chart AND download buttons actually
                                work end to end (real browser downloads, not
                                just "menu item exists" - this is what caught
                                the _oraChartId/hardcoded-div-id mismatch that
                                silently broke their downloads), that the
                                ULM regression chart's legend does not overlap
                                its xaxis title on the pane's FIRST reveal (a
                                real, user-reported bug: this chart is built
                                eagerly while its pane is still display:none,
                                which drew the legend ~20px too high; the old
                                fix, Plotly.Plots.resize() on reveal, rescales
                                the plot but does not recompute a
                                miscalculated legend position - confirmed by
                                also trying relayout({}) and Plotly.react(),
                                neither of which moved it either - only a
                                fresh Plotly.newPlot(), now stashed per-chart
                                as div._omiscerRedraw and used by activate(),
                                fixes it; checked here BEFORE the TF-switch
                                checks below, which redraw the chart from
                                scratch and always incidentally masked the
                                bug), and that no unexplained console errors
                                occurred.
playwright_genomics_results_cards.py  13 checks. The FIRST browser coverage
                                either Genomics 'Network Propagation / Module
                                Detection' or 'Gene Function Annotation' has
                                ever had (no genomics Playwright file existed
                                before this one). Same real, user-reported
                                box-wrapper bug as playwright_regulon.py
                                above, fixed the same way, for these two
                                screens - drives the REAL upload -> resolve
                                -> configure -> run flow for each (no
                                synthetic STATE injection, which would risk
                                silently drifting from the real payload shape
                                and giving false confidence).
playwright_transcriptomics_network.py  41 checks. The FIRST real-browser
                                coverage of the Transcriptomics Network /
                                Module Detection tile - the fourth
                                Transcriptomics analysis, replacing the
                                disabled "Enrichment (coming soon)"
                                placeholder. Drives the real flow (tile ->
                                organism selector scoped to STRING, not the
                                whole UniProt taxonomy -> resolve -> config ->
                                run -> results) at 500 permutations and
                                asserts actual CONTENT, not just that panes
                                exist: transcriptomics-shaped hit columns and
                                the null-distribution figure rendering. Under
                                the corrected, production-scale-calibrated
                                statistics, this specific real run legitimately
                                produces 0 significant proteins - so hit/
                                recruited/caveat/module/network-figure
                                assertions are REPORT-ONLY or gated on the
                                real observed count (stats["hits"]>0,
                                stats["seedHits"]>0, stats["recruitedAll
                                Family"]>0) rather than assuming a positive
                                result, keeping the suite honest either way
                                rather than loosened. Also covers three
                                evaluation fixes verified on this SAME real
                                run: (3a) a genuinely underflowed FDR renders
                                "< 1e-300" on screen, never the pre-fix
                                "0.000e+0", while the CSV download keeps full
                                precision (asserted unconditionally, not
                                data-dependent); (3b) the "recruited only" Hit
                                table toggle's default state matches whether
                                recruited genes actually exist
                                (toggle_state == stats["recruited"]>0), with a
                                full branch for the zero-recruited case
                                (verifies default OFF, full table shown, and
                                that checking it correctly empties the table
                                to 0 rows) and, symmetrically, would verify
                                the positive-path filtering if a run produces
                                recruited genes; (3c) the on-screen caveat's
                                own parsed numbers match the payload's
                                n_hits_seed/n_hits/n_recruited_all_family
                                exactly, each check independently gated on
                                whether that specific caveat is expected to
                                appear at all. PLUS (multi-database): the
                                config screen's "Label modules using" row is
                                now checkboxes (GO + every MSigDB collection
                                pill, since renderMsigdbCollectionSelector
                                reuses the same ora-db-pill markup) rather
                                than a single dropdown, GO is checked by
                                default, MSigDB starts unchecked, and the
                                draft's own enrichmentDatabases is confirmed
                                to be an ARRAY (["GO"]), not the old single
                                string.
playwright_proteomics_ora.py    57 checks. The FIRST real-browser coverage of
                                the Proteomics module (Milestone 1: the
                                foundation plus Over-Representation
                                Analysis). Drives the real flow end to end -
                                tile -> sub-tile gating (ORA and, since
                                Milestone 2, Quantitative Enrichment both
                                enabled; PPI Network / PTM confirmed enabled
                                too once those milestones landed - see
                                playwright_proteomics_network.py and
                                playwright_proteomics_ptm.py for their own
                                dedicated coverage) -> organism -> upload -> REAL
                                UniProt ID Mapping resolve of all 3500
                                example accessions -> configure -> run
                                against real InterPro/Complex Portal/UniProt
                                fetches -> results. Asserts biology, not just
                                that panes exist: Complex Portal's top hit
                                must be one of the six real complexes
                                actually seeded into the bundled example.
                                Also covers the messy MaxQuant example's
                                contaminant count reaching the SCREEN (not
                                just the payload); DISEASES being present,
                                enabled for human and correctly greyed out
                                for mouse; the custom-DB example zip's exact
                                contents downloaded through the real
                                endpoint, then re-uploaded and parsed; a
                                symbol-based custom DB being REJECTED rather
                                than silently matching nothing; and a
                                sensitivity check that runs the real
                                _oraGeneSet against a real result row to
                                prove the renamed "Overlapping Proteins"
                                column (formatted "ACCESSION (SYMBOL)") is
                                still parsed into bare accessions - without
                                it the protein-concept network's log2FC
                                colouring would silently stop resolving.
                                Also drives KEGG through the real UI twice:
                                a clean human run (KEGG pill present, no
                                organism picker, real terms tested), then a
                                fresh upload/resolve for the real Bacillus
                                subtilis ambiguity, where the SAME
                                disambiguation-picker widget Transcriptomics
                                uses (generalized with moduleKey="proteomics")
                                appears, offers the real 14-strain candidate
                                list, and confirming persists into
                                Proteomics' OWN state key
                                (STATE.proteomics._keggOrganismConfirmed, not
                                a Transcriptomics one) - proving the
                                generalization genuinely works for a second
                                module, not just a Transcriptomics-shaped stub.
playwright_proteomics_quant.py  56 checks. Milestone 2 - the first real-
                                browser coverage of Quantitative Enrichment
                                (1D annotation enrichment / competitive
                                Wilcoxon, reusing Genomics' Gene-set Burden
                                Test engine via alternative="two-sided").
                                Drives tile -> organism -> upload -> REAL
                                UniProt ID Mapping resolve -> continue ->
                                configure (asserts NO padj/log2FC cutoff
                                fields on this screen, and that the ranking-
                                metric selector defaults to log2fc) -> run
                                against InterPro/Complex Portal/UniProt
                                keywords -> results, asserting result keys
                                carry NO direction suffix (unlike ORA's
                                "(all)/(up)/(down)") and no ORA-style
                                significance-count field leaked into the
                                summary. The load-bearing sensitivity check:
                                _burdenGeneSet - which read ONLY row['Genes']
                                before this milestone - is run against a REAL
                                Proteomics row (column renamed to "Overlapping
                                Proteins") and must return a NON-EMPTY,
                                bare-accession set; before the fix it would
                                have silently returned empty, breaking the
                                enrichment map and protein-concept network's
                                edges with no error. Also selects whichever
                                selected database actually came back
                                significant this run (Complex Portal,
                                live-verified 2026-08-25 after fixing the
                                bundled example's per-protein-random-
                                direction generation bug - see
                                generate_proteomics_examples.py) and
                                exercises all 5 Burden-shaped
                                chart builders (dot, bar, volcano, enrichment
                                map, protein-concept network) for real
                                rendered Plotly figures, plus the CSV download
                                endpoint. Also closes a real gap found while
                                auditing the module for completeness: MSigDB
                                was selectable here (same as ORA/GSEA) but had
                                NO cross-collection summary figure - ORA and
                                GSEA both have one, Quant did not. Runs a real
                                MSigDB: Oncogenic pass, opens the new summary
                                figure and confirms it renders with REAL data
                                (proving _msigdbFamiliesFromQuantTables' own
                                column names - "FDR", "Effect (mean
                                difference)" - actually match what the
                                Wilcoxon engine's blocks carry, not just that
                                the function is defined), then changes the
                                significance threshold with the summary still
                                open and confirms it survives - the one
                                behaviour easy to miss when wiring
                                refreshMsigdbSummary (declaring it without
                                actually calling it from the alpha input's
                                change handler). Also drives a clean-organism KEGG
                                run through this screen's own OWN checkbox and
                                widget instance (the full strain-picker flow
                                is already covered in
                                playwright_proteomics_ora.py, using the same
                                generalized widget - no need to duplicate it
                                here), confirming a real single "KEGG" family
                                with no direction suffix comes back with
                                tested terms.
test_organism_allowlist_mirror.py  7 checks. Guards a drift risk found while
                                auditing organism generality: several
                                per-database organism allowlists are
                                DUPLICATED between Python (the source of
                                truth) and JS (the UI mirror) - if Python
                                ever gains a taxon (e.g. TRRUST starts
                                publishing rat), the JS mirror will not pick
                                it up automatically, and nothing today would
                                say so. Compares the two live against each
                                other - TRANSCRIPTOMICS_REGULON_TAXA vs
                                OMNIPATH_TAXA/TRRUST_URL_BY_TAXON,
                                _regulonDbAvailableForTaxon vs each regulon
                                source's actual per-source gating (exhaustive
                                over every database x a probe set including
                                rat and Arabidopsis, which Python correctly
                                rejects for all four), _msigdbAvailableForTaxon
                                vs MSIGDB_TAXON_SUFFIX, and
                                _msigdbCollectionAvailableForTaxon vs
                                msigdb_collection_available for every
                                collection - rather than against a hardcoded
                                literal, so a Python-side change is caught
                                automatically. A sensitivity check confirms a
                                deliberately-broken "always available" gate
                                genuinely produces mismatches against the
                                real Python source of truth.
test_gsea_fdr.py               24 checks. Correctness of GSEA's two FDR/
                                q-value columns (Transcriptomics_module.py,
                                also shared by the Regulon activity module).
                                Written because _pooled_fdr's own docstring
                                claimed its monotonicity "is verified
                                explicitly in this module's unit tests" and no
                                such test existed anywhere in this directory -
                                an audit question ("is the FDR q-value
                                calculated correctly?") surfaced the gap.
                                _pooled_fdr (the GSEA-paper pooled-permutation
                                FDR, not a plain BH correction on nominal p)
                                is cross-checked against an INDEPENDENT,
                                literal brute-force reimplementation of the
                                paper's formula (plain loops, no shared
                                helpers with the vectorized code under test) -
                                exact match (0.0 max diff) across 5 random
                                seeds. Also checks monotonicity within each
                                NES-sign family directly (q must not increase
                                as |NES| gets more extreme), with a
                                SENSITIVITY check that the compared sequence
                                genuinely varies first (a constant sequence
                                would pass "monotone" vacuously either
                                direction). Confirmed non-vacuous by actually
                                injecting the exact "backward running-min"
                                bug the real function's own docstring warns
                                about (wrong-direction monotonicity) via
                                runpy in-process monkeypatching - 9 of the 24
                                checks go red, including 5 of the 5 brute-
                                force cross-checks. Edge cases: an entirely
                                one-signed scenario (no negative-NES terms at
                                all), and a null pool with no same-sign values
                                for one term's sign (denom_null_total=0) falls
                                back to q=1.0 rather than NaN or an unhandled
                                division. benjamini_hochberg (the sibling
                                "FDR (BH)" column, classic BH on nominal p) is
                                cross-checked the same way against an
                                independent reference across 3 seeds, plus the
                                two BH invariants (q non-decreasing in
                                p-rank; q never below its own p) and the
                                empty-input edge case.
test_gsea_null_histogram_range.py  4 checks. A real, user-reported bug in
                                buildGseaNullHistogram (static/app.js): when
                                the observed NES is more extreme than the
                                ENTIRE permutation null distribution's own
                                range - exactly what a strongly significant
                                result looks like - Plotly's default
                                autorange stretched the x-axis just far
                                enough to include the observed-NES line, with
                                near-zero padding, pinning it onto the axis
                                border where it was visually indistinguishable
                                from it (confirmed against a real screenshot:
                                NES=-2.40, p=6.8e-3, no line visible at all).
                                Fixed with an explicit xaxis.range padded 6%
                                past whichever is more extreme - the null
                                samples, the observed NES, or 0. A sensitivity
                                check reverted the real fix and confirmed the
                                exact same test goes red against the
                                unfixed code.
test_msigdb_summary_figure.py   15 checks. The MSigDB cross-collection
                                summary figure, driven in a real browser
                                against SYNTHETIC tables (no network). Exists
                                because a post-implementation audit found
                                three defects the "it rendered with data"
                                check in playwright_msigdb.py could not see:
                                the bar colour applied a signed rule to ORA's
                                UNSIGNED fold enrichment (making the colour
                                channel vacuous), a non-numeric FDR in the
                                first row pinned the family's "best hit"
                                forever via NaN comparison, and a null FDR
                                coerced to 0 (i.e. maximally significant).
                                Two sensitivity checks re-run the assertions
                                against the pre-fix logic to prove they are
                                not vacuous. Also covers a real, user-reported
                                display issue: scipy's hypergeom.sf() can
                                genuinely underflow a p-value (and hence its
                                BH-FDR) to a literal 0.0 for an astronomically
                                significant term, and the old ".toExponential(2)"
                                formatting rendered that as "0.00e+0" -
                                indistinguishable from the null-FDR-coerced-
                                to-0 bug above, even though it is a real,
                                correct result. Now shown as "< 1e-300"; a
                                sensitivity check reverted the real fix and
                                confirmed the check goes red with the exact
                                misleading string a user actually saw.
playwright_transcriptomics_kegg.py
                                45 checks. Real-browser flow for KEGG in
                                Transcriptomics ORA and GSEA, and the only
                                machine check on static/app.js in an
                                environment with no Node (any console or page
                                error fails the run). Covers what neither
                                file above can see: the KEGG pill exists on
                                BOTH config screens, checking it triggers a
                                full re-render so its organism caption
                                actually appears (a real bug this test caught
                                - the toggle only re-rendered the run button),
                                and - the point of the whole design - NO
                                organism search box is rendered for KEGG,
                                since a second organism control here is
                                precisely the divergence bug of commit
                                27ce033. Then runs both analyses for real:
                                three "KEGG (direction)" families whose term
                                IDs all start with "hsa" (proving the DERIVED
                                organism was really the one used), a slugged
                                result card, the same four result-nav panes
                                every other database gets, and a figure that
                                actually draws when its lazy "Show ..." toggle
                                is clicked.
                                Also drives the strain-disambiguation picker
                                through the real UI against the real Bacillus
                                subtilis ambiguity (not a stub): checking KEGG
                                for that organism makes the picker appear with
                                the real 14-strain candidate list by name,
                                confirm stays disabled until one is chosen,
                                confirming replaces the picker with a
                                "Confirmed manually" caption and persists the
                                choice into client state, a real run then
                                carries the manual-confirmation caveat in its
                                result table (proving the confirmed organism
                                genuinely reached the run, not just the UI),
                                and switching back to a clean organism clears
                                the confirmation from client state.

playwright_msigdb.py            68 checks. Same real-browser approach for
                                MSigDB in ORA and GSEA: the collection
                                sub-selector only appears once "MSigDB" is
                                checked, mouse organism greys out exactly
                                PID and Oncogenic (confirmed live to be the
                                only two collections genuinely missing a
                                mouse file - WikiPathways does NOT belong in
                                that list, an earlier draft of this test got
                                that wrong), a real run (Oncogenic, kept
                                small for speed) produces the new
                                "MSigDB: <collection> (direction)" result
                                keys, the family-count caption appears, and
                                the MSigDB collection summary figure renders
                                with real data - in both ORA and GSEA. A
                                sensitivity check confirmed the mouse-gating
                                assertion goes red if
                                _msigdbCollectionAvailableForTaxon is broken
                                back into a blanket "true". Also covers the
                                two audit fixes that need the real UI: the
                                summary figure rebuilding when the FDR
                                threshold changes (it used to cache the first
                                alpha forever, showing bar lengths and an axis
                                title contradicting the table below it), and
                                the organism-consistency fix: switching
                                organism after a resolve clears the resolved
                                table and every result, shows a warning naming
                                the new organism, and - caught by this test
                                itself - does NOT leave a stale, now-wrong
                                warning on screen if you switch again before
                                resolving (an earlier version of the fix only
                                set the notice, never cleared it, so a second
                                switch back to the original organism kept
                                claiming it had just changed to the first
                                one). And result-table CSV downloads carrying
                                the result key in the filename (every family
                                previously saved as one fixed name, so a
                                multi-collection run produced a pile of
                                indistinguishable "... (1).csv" files).

playwright_metabolomics.py      423 checks. Real-browser coverage of the five
                                Metabolomics tiles, written because the first
                                Metabolomics UI was never opened in a browser
                                and was rejected on sight (statistics columns
                                metabolomics does not use, every database on
                                one screen, a look unlike the rest of the app).
                                Asserts: five separate tiles; no log2FC/pvalue/
                                padj input on the list tiles; house-style
                                structure - flat module screen with the shared
                                file widget and inline example link (no
                                numbered card-per-step wizard), exactly ONE
                                result <select> and ONE result card; every
                                chart control present AND redrawing through
                                each option; and chart hover/axis text saying
                                "Compound", never "Gene" (the builders were
                                written for genes). Also the Quantitative
                                tile's layout/group override round trip, the
                                derived statistics shown and downloadable, and
                                burden columns renamed for compounds.
                                Also: checking KEGG on a pathway tile shows the
                                KEGG organism box (the ONLY organism search
                                box on the setup screen), which reports the
                                derived organism by default; typing
                                'saccharomyces' searches real KEGG and lists
                                'code | name', picking validates it and can be
                                cleared, and on the Quantitative tile an
                                explicit KEGG organism beside a human organism
                                reaches the run (sce maps, caveat says chosen
                                explicitly); unchecking KEGG removes the box;
                                the CLASS tiles offer NO enrichment map
                                (asserted absent - class terms are disjoint,
                                so the map could never have an edge) but keep
                                the dot plot, and DO offer a class-composition
                                pie:
                                its slices sum to the number of distinct
                                matched compounds, one colour picker per real
                                class, a picked colour recolours the drawn
                                slice and survives a threshold change, the
                                level and "Classes shown" controls redraw
                                (capped, with nothing lost), every slice has a
                                distinct colour. The pie's percentage labels:
                                Plotly's own slice text is off (it drops small
                                slices) and every slice gets a label as a
                                draggable annotation; SENSITIVITY - each
                                leader line must end inside ITS OWN slice just
                                inside the rim (checked from the angles, at
                                100% and 200% zoom; a deliberate mutant that
                                mirrored the tips turned 6 checks red);
                                Hide/Show, a real mouse drag moves a label
                                while its line stays attached, the position
                                survives a redraw, Reset restores it, and Zoom
                                200% doubles the figure. Slice borders (white /
                                black / none) and "Label slices of at least"
                                (a threshold chosen per panel so it genuinely
                                separates that panel's slices). A synthetic
                                60-CLASS pie sizes itself: every slice
                                labelled, no two labels closer than a line, none
                                clipped, the whole 60-entry legend fits without
                                scrolling, Zoom 50%/200% keep every label.
                                REGRESSION (found in a 1280px window): long
                                class names made Plotly grow the right margin
                                for the legend, shrinking the plot area and
                                sliding the pie out from under its leader tips;
                                asserted at 700/900/1280px that the plot area
                                keeps the size asked for and every tip sits on
                                the rim of the pie Plotly really drew. A mutant
                                that removed the legend-width reservation
                                collapsed the plot area to 76px and turned
                                those checks red. No pie on any pathway panel,
                                including the KEGG panel of the mixed-family
                                Quantitative tile.
                                Pathway Topology Analysis: the input-mode
                                choice (defaults to "list"), the node-
                                importance selector and KEGG organism box on
                                setup, NO database picker (there is only
                                ever one), and on results NEITHER the usual
                                "Compound set overview" pane NOR a class-
                                composition pane - instead a "Pathway impact"
                                pane and a "Pathway graph" pane. Also
                                concentration-table mode: switching the input
                                select updates topologyInput, the same
                                group-contrast fields Quantitative Enrichment
                                uses appear, the chosen measure and Wilcoxon/
                                burden column shape reach the run, and
                                topology's own statistics preview is shown.
                                Pathway impact plot: a categorical dot plot
                                (buildDotPlot's own convention) with every
                                pathway's real NAME on the y-axis, not
                                floating text over a continuous axis -
                                checked that the trace stays mode "markers"
                                (no on-plot text of its own) and every y
                                value is a real name; EVERY tested pathway is
                                drawn by default (point count == row count,
                                not just the significant ones), and "Terms
                                shown" narrows it (asked 10, got exactly 10).
                                Pathway graph: node markers carry a real
                                COMPOUND name (not a bare KEGG id - checked by
                                regex), but as a Plotly ANNOTATION anchored to
                                each node, not the marker trace's own text -
                                same mechanism as the class-composition pie's
                                percentage labels (buildClassCompositionPie),
                                reused here for the same reason: a SENSITIVITY
                                check drags one node's label with a REAL mouse
                                gesture and confirms the annotation's anchor
                                (x, y - the node's own KGML position) does NOT
                                move while its tail (ax, ay - the label) does,
                                the moved position survives a Hide/Show
                                redraw (stored on the run, like the pie's
                                dragged labels), and "Reset label positions"
                                puts it back; "Compound labels" -> Hide drops
                                every annotation (0 left) rather than just
                                changing trace mode. Annotation/plot queries
                                are taken off the LAST .js-plotly-plot in the
                                DOM, since the impact pane's chart div is
                                hidden, not removed, once shown; the caption
                                explaining isolated (unconnected) nodes as
                                genuine KEGG data, and the one explaining
                                nodes/labels can be dragged, are both present.
                                Impact plot REGRESSION (found from a user
                                screenshot): a hardcoded x-axis range clipped
                                the biggest marker at either end until
                                "Autoscale" was clicked by hand - checked that
                                xaxis.autorange is genuinely true on first
                                draw (gd.layout.xaxis.range itself is NOT a
                                reliable signal: Plotly writes the computed
                                range back into it once autoranging runs, so
                                it is never undefined post-render). Node drag:
                                _makeNetworkDraggable (already used by the
                                enrichment map / gene-concept network charts)
                                was wired onto the pathway graph's own nodes
                                too, with an onMove callback that keeps a
                                dragged node's OWN label annotation anchored
                                to its new position - checked with a real
                                mouse gesture (node position computed the same
                                way the helper's own hit-test does, since a
                                Plotly marker has no selectable DOM element)
                                that the node moves AND its annotation's x/y
                                anchor moves with it, to the same new spot.
                                SENSITIVITY (a real regression hit while
                                building this): node-dragging's hit-test
                                tolerance (20px) reaches a label sitting only
                                16px from its own node by default, so it was
                                silently stealing every click meant for
                                Plotly's own annotation-tail label drag,
                                breaking the label-drag checks above the
                                moment node-dragging was added - fixed by
                                having _makeNetworkDraggable's mousedown
                                handler step aside (evt.target.closest
                                ('.annotation')) for a click that lands on an
                                annotation, a no-op for every OTHER chart that
                                uses the same shared helper since none of them
                                carry annotations.
                                Needs network.

playwright_multiomics_mapping.py  105 checks. Real-browser test of Cross-omics Pathway
                                Mapping on real data (port 8009, real network): the
                                second Multi-omics box is live and shows only its name;
                                four omics -> upload -> resolve -> 'Continue to pathway
                                mapping' -> the KEGG picker (340 pathways, capped at
                                150 rows, search by name and by ID) -> KEGG hsa00020 ->
                                the four panes -> a pathway the expression tables DO hit
                                (direction triangles whose shape agrees with the
                                direction in their own hover) -> Reactome R-HSA-71403 and
                                the diagram-less Glycolysis (the overview says it was
                                drawn inside Glucose metabolism) -> back, and re-opening
                                the box reuses the resolved uploads (the Joint box still
                                starts clean) -> 'Map this pathway' from a Joint result.
                                THE EVIDENCE CHIPS NEVER COVER EACH OTHER: the pixel
                                distance between every pair of chips on one box is
                                measured at first draw, zoomed in and zoomed back out
                                (11.0px against a 9px marker), on a map with boxes
                                carrying three or more chips - without such a box the
                                check would be vacuous. Downloads are real CSVs.
                                SECOND ROUND, in the real browser: hsa01100 draws with
                                no pile in a corner (0 of 1,342 handles within 4px of
                                the top-left-most one), grey lines for the genes nobody
                                hit and coloured ones for the hits, only the hits
                                labelled by default (the count equals the number of hit
                                boxes); hsa04210 draws its 135 relations, the dark ones
                                are exactly the relations with both ends hit, 4 group
                                outlines, an 'Interaction links' pane and CSV, and
                                every hit box has its chips; Reactome's Disease shows NO chart but the
                                sub-pathway table, and pressing Map down the nesting
                                reaches a real diagram.
                                THIRD ROUND: hsa01100 is drawn at its own size (about
                                5,000 x 3,300px, in a 1,050px frame that scrolls), has a
                                colour picker per pathway class, every hit gene is a
                                thick full-strength line with its chips, choosing red for
                                Lipid metabolism recolours exactly its lines and
                                'KEGG colours' restores them, and the label overlap is
                                MEASURED: 154 of 204 hit labels sit on another label when
                                the map is squeezed to the page, 21 of 204 at 100%. The
                                zoom tools are checked on this real 5,000px map too.
playwright_multiomics_mapping_chart.py  114 checks. The same chart on SYNTHETIC data, in
                                seconds (the server only serves app.js): a hand-built
                                result is put into STATE and the real render path draws
                                it, so every expected value is known in advance, at two
                                window widths. Shapes (triangle up / down, diamond for
                                disagreeing members, filled circle for a list hit, a
                                HOLLOW circle in its own trace for measured-not-
                                significant, no chip when absent), FIXED slots (an
                                absent omics leaves its slot empty instead of shifting
                                the others), chips below their box, no label on a box
                                without evidence, a darker border on a hit, chip
                                spacing >= a marker at first draw / zoomed in 5x /
                                zoomed out / far out, KEGG substrate-to-product lines
                                vs Reactome's own segments with dotted catalyst links,
                                links to other pathways as plain text, empty-state
                                messages. 19 frontend mutants (stacked chips, no
                                re-placement, no relayout listener, swapped triangles,
                                filled instead of hollow, chips above the box, ...):
                                17 caught; the other 2 change no behaviour (a list omics
                                never has a direction, so 'draw a list hit as a
                                triangle' is the same code path, and a link node never
                                reaches the box traces). app.js restored byte-identical.
                                SECOND ROUND: the relation styles (arrow, bar, plain,
                                dashed arrow, dash-dot; dark only when BOTH ends are hit;
                                an arrow stops short of its target; '+p'), group
                                outlines, label modes (the default on a map of over 300
                                boxes is hits only), the overview screen, and a teardown
                                regression check (not mutation-sensitive: see its
                                comment). THIRD ROUND (dragging was taken out at the
                                user's request: it was incomplete and hard to use):
                                line genes coloured by pathway CLASS (hit = full
                                strength and thick, measured = firmer, untouched =
                                faint, hits on top, drawn along the polyline with y
                                flipped, chips at a hit line's middle, no straight
                                reaction lines, pathway labels in their class colour); a
                                colour picker per class that recolours exactly that class
                                (and its labels) and 'KEGG colours' that restores it; a
                                wide map defaults to 100% in a scrolling frame, the Size
                                menu scales the chart exactly (50% / 200%), 'fit' is the
                                page's width, the wheel scrolls at a fixed size, chips
                                stay a marker apart at every size, 'Download map' opens
                                a menu you can reach on a 2,160px chart; and dragging
                                moves nothing. 25 more mutants over the round, all
                                caught; app.js restored byte-identical. TOOLS: at a fixed
                                size the chart's own toolbar is off-screen, so Zoom tool,
                                Pan tool, Zoom in / out, Autoscale and Reset axes are
                                repeated above the chart (hidden at fit-to-page): Zoom in
                                halves the visible range about its centre and re-places
                                the chips, Zoom out doubles it, Reset axes returns EXACTLY
                                to the opening view (Plotly widens the x range to the
                                canvas's shape, so that, not the raw data box, is the home
                                view - my first assertion compared against the wrong one),
                                Autoscale after zooming in contains every box again, the
                                two tools switch the drag mode. 8 more mutants, all caught.
playwright_multiomics.py        127 checks. Real-browser test of Multi-omics Joint Pathway
                                Enrichment (static/app.js): the formerly INERT
                                'Multi-Omics' button -> choose omics -> one upload
                                each -> resolve -> setup -> run -> results -> a
                                second, four-omics run. SENSITIVITY: clicking the
                                button opens the screen (it used to do nothing);
                                no upload widget appears until TWO omics are
                                chosen, and only for the chosen ones; Stouffer is
                                NOT in the method select; DE-cutoff fields appear
                                only when a DE omics is chosen; the joint dot plot
                                is a categorical dot plot with pathway NAMES on a
                                reversed y-axis (Plotly stores 'reversed' as
                                autorange true + a high-to-low range, so the check
                                reads the range); the evidence plot has one trace
                                per omics plus a combined diamond, in each omics'
                                own colour, and the metabolomics 'Gene' ->
                                'Compound' chart-noun swap does not mangle the
                                omics names. REGRESSION found by this test: a
                                member list of 82 characters sits under
                                buildTable's 100-character clamp but still wrapped
                                to ten lines in the ~100px column of an 18-column
                                table, inflating its row to 195px - every member
                                cell is now clamped to three lines (full text in
                                the tooltip and the CSV). Test-harness lesson:
                                Playwright fill + a synthetic change leaves the
                                browser's own change pending, and it re-renders the
                                panel under the NEXT click - commit with Tab. Every
                                chart control redraws with no page error; the
                                downloaded CSV has every combined row; the active
                                pane survives a database switch. Mutation-tested
                                for the row clamp only (removing it turns the
                                195px-row check red); the other rules are
                                asserted but were not individually mutated.
                                TILE SCREEN: Multi-Omics opens six boxes (the
                                planned methods); only Joint Pathway Enrichment is
                                clickable, the other five are disabled and say
                                'coming soon', and Back goes one level at a time.
                                LANES: the evidence plot splits each pathway row
                                into one lane per omics plus the diamond, because
                                two omics at the same p-value (an omics with no hit
                                sits at exactly 0) used to stack and hide one
                                marker; checked in pixels on the four-omics data -
                                coincident markers stay >= 9px apart (12.9px
                                measured; 0px with the lanes collapsed). Needs
                                network. A full run takes ~12 minutes (the
                                four-omics step resolves ~6,500 identifiers).

These six (playwright_regulon.py, playwright_genomics_results_cards.py,
playwright_msigdb.py, playwright_transcriptomics_network.py,
playwright_proteomics_ora.py, playwright_proteomics_quant.py) make real
network calls (mygene.info, TRRUST, ChEA3, MSigDB, UniProt, STRING,
GO/QuickGO, InterPro, Complex Portal, DISEASES) and are not mocked, so they
need internet and take a few minutes. test_msigdb_summary_figure.py,
test_gsea_null_histogram_range.py, and test_organism_allowlist_mirror.py
are the exceptions - all three need the server only to load app.js, and
make no external calls.


---------------------------------------------------------------------------
RUNNING EVERYTHING
---------------------------------------------------------------------------
There is no runner script - each file is independent by design, so a failure
in one never hides another. From the repo root, with the server up:

    Get-ChildItem tests\*.py | ForEach-Object { python $_.FullName }

Expected totals as of the last full run, in file-listing order: 117, 22, 16,
4, 15, 5, 16, 10, 24, 58, 44, 49, 52, 8, 23, 19, 20, 31, 23, 31, 49, 21, 15,
12, 12, 22, 22, 16, 12, 49, 39, 90, 13, 37, 40, 41, 7, 4, 15, 68 - all passing
(1171 checks across 40 files). The Metabolomics files listed above were added afterwards
and are not in that tally.

MULTI-OMICS CORRELATION-BASED INTEGRATION (info/README-MultiOmics-Correlation.txt).
Run with the server on port 8009 only for the live and browser files.

test_multiomics_correlation.py  105 checks. The engine (Multiomics_correlation), offline.
                                Every statistic against INDEPENDENT pure-python code:
                                average ranks, Pearson and Spearman with ties and
                                missing values (each pair with its own n; Spearman
                                RE-RANKED over each pair's shared samples, as R and
                                pandas do - a check proves the test can tell this
                                from ranking each feature on its own), a
                                hand-written incomplete-beta p-value (no scipy t), and
                                Benjamini-Hochberg against Genomics_module's own.
                                NULL CALIBRATION on independent random matrices: about
                                5% of p-values below 0.05, a Kolmogorov-Smirnov bound
                                for uniformity, nothing surviving BH (n = 20 and 12).
                                Planted pairs recovered; sample alignment (whitespace,
                                the case option both ways, the ones left out named, the
                                'other orientation' hint); orientation detection both
                                ways with a Group column, and annotation columns
                                (gene_biotype, a compound Class, 'Compound group')
                                that must NOT flip it; the 20% missing limit; log2;
                                duplicate-gene averaging; the clustering order against a
                                naive average-linkage; the layered layout against the
                                input order. No np.corrcoef / np.cov / scipy spearmanr:
                                they hard-crash this environment, so the engine uses
                                elementwise numpy only.
test_multiomics_correlation_session.py  73 checks. The routes and jobs in API.py,
                                offline with a stubbed resolver: the order of the
                                steps, Genomics refused, fewer than two omics refused,
                                a re-upload / a re-choice / a re-resolve drops what
                                depends on it, a failed run leaves no half-built result,
                                the example files end to end, the scatter endpoint (its
                                axis labels once overwrote the value arrays - a
                                regression check; a protein axis names gene and
                                accession), a repeated id 'CS (2)' resolved through
                                'CS', downloads, strict JSON, isolation
                                from Joint and Mapping.
live_e2e_multiomics_correlation.py  32 checks. Real HTTP end-to-end (port 8009, real
                                resolvers). Rebuilds from the raw CSV files, with no
                                Omiscera code: every gene's rho / n / p / FDR in the
                                concordance table, and the COMPLETE transcript x
                                metabolite block (12,558 pairs, BH): the reported set
                                must be identical to the app's; the same for the
                                protein x metabolite block, whose proteins have
                                missing values (every pair ranked over its own
                                samples, one pair at a time). Planted structure
                                (CS-citrate, IDH3A-isocitrate, MDH2-malate; the
                                discordant genes negative); downloads equal the tables.
playwright_multiomics_correlation.py  78 checks. Real browser (port 8009): the third
                                box is live; Genomics cannot be chosen; the orientation
                                menu re-reads the same file (forcing the wrong layout
                                visibly changes the sample count); the alignment card;
                                the concordance histogram, table and a scatter whose
                                rho equals the table's; pairs and their scatter; the
                                heatmap (symmetric colour scale, a dot on exactly the
                                reported pairs, a fixed cell size measured on a 40 x 40
                                synthetic heat); the network (every node drawn, no two
                                nodes of a column overlap, middle-column labels above
                                their nodes); the figure tools row on every chart;
                                the Setup block has no overlapping controls (measured
                                boxes); the colour-gradient menu (6 diverging scales +
                                swap) reaches the heatmap, the histogram and the
                                network, and the caption never names a colour;
                                changing the case option drops a stale match; a
                                2,500-row pair list draws 2,000 rows with a note
                                (the number format follows the browser's locale);
                                downloads; no page error. Set CORR_SHOT_DIR to save a
                                screenshot of each pane to look at.
generate_correlation_examples.py  Regenerates Data3/example_corr_*.csv (synthetic, seeded,
                                planted structure; needs the network once for the
                                UniProt accessions). Not a test.

MULTI-OMICS NETWORK-BASED INTEGRATION (info/README-MultiOmics-Network.txt).
Run with the server on port 8009 only for the live and browser files.

test_multiomics_network.py      79 checks. The engine (Multiomics_network), offline, on a
                                3-community STRING-like fixture. The equal-weight rule against
                                hand computations (one omics = exactly 1.0; 3,000 transcript seeds
                                and 4 metabolite seeds carry equal total weight; a protein hit
                                by two omics outweighs one hit by one). The metabolite -> enzyme
                                join against a brute-force join on random data at several limits,
                                the InChIKey-skeleton bridge, the hub guard at its boundary
                                (60 reactions kept, 61 left out) and the family-reaction guard at
                                its boundary (30 proteins kept, 31 skipped). THE KEY CHECK:
                                identical seeds from two omics reproduce the existing Proteomics
                                network tile's propagated score, Z, p, FDR and modules exactly
                                (difference 0.0). Seeds named by omics, per-omics module counts,
                                the DE columns travelling with a protein. The ChEBI skeleton
                                cache remembers 'no structure' but never a failed lookup. The
                                pathway network against a brute-force Jaccard (threshold boundary
                                inclusive, a compound named like a gene kept apart, components,
                                top-N) and a layout with no two nodes closer than the gap.
test_multiomics_network_session.py  60 checks. The routes and the job in API.py, offline with a
                                synthetic STRING graph and stubs: the order of the steps, an
                                organism STRING does not have (refused BY NAME), every parameter's
                                range, nothing half-valid stored, the strict-JSON payload, four
                                omics in one run, the metabolite table, both networks, the
                                downloads, isolation from the single-omics tiles and from Joint,
                                Mapping and Correlation, clearing on re-upload / re-resolve.
live_e2e_multiomics_network.py  51 checks. Real HTTP end-to-end (port 8009, real STRING human network at
                                confidence 700, real Reactome files, 500 permutations). Rebuilds
                                from the raw files with no Omiscera code: each omics' seeds (STRING's
                                own names and aliases; proteomics' padj by Benjamini-Hochberg), the
                                metabolite -> enzyme join (Reactome's two files streamed again),
                                the combined weights, and THE RANDOM WALK itself on the rebuilt
                                graph - every protein's propagated score must match. Also the
                                module seed counts and the pathway edges.
playwright_multiomics_network.py  41 checks. Real browser (port 8009, real data): the fourth box is live;
                                the shared upload screen leads to the network setup; the setup
                                shows the metabolite limits only because metabolomics is chosen;
                                the overview, the recruited-only hit table, the single-omics
                                tiles' own figures drawn from this payload, the module
                                composition chart, the metabolite table, and the pathway
                                network at a fixed readable size in a scroll frame: every node
                                drawn, no two nodes too close, NO TWO LABELS OVERLAP (measured
                                on the page), one chip series per omics, the figure tools row.
                                Set NET_SHOT_DIR to save a screenshot of each pane.
playwright_multiomics_network_chart.py  58 checks. The same results screen on a SYNTHETIC payload
                                that HAS what the real example lacks (seconds, no run): recruited
                                proteins and their default view, modules with seeds from several
                                omics, the composition chart (stacked, bars equal the table), the
                                single-omics tiles' figures, module enrichment, the metabolite
                                hub row, and a crowded pathway network (35 long Reactome names
                                laid out by the engine): positions equal the engine's, NO TWO
                                LABELS OVERLAP, labels show no visible '<br>', chips only where
                                an omics hit and under their node AS DRAWN, near-identical
                                pathways merged and said so in the hover. THE FIGURE TOOLS: a
                                visible tools row (zoom, pan, zoom in/out, autoscale, reset axes)
                                that works on the chart (zoom in halves the view, reset returns
                                to the opening view), Size (50-200% and fit to page; chips keep
                                their place at every size), name density, and the Download
                                figure menu opening under its button. Every rule goes red
                                when its line of app.js is broken.

A NOTE ON WRITING NEW CHECKS HERE. Prefer assertions that can actually fail
for the reason you care about. A concrete example lives in
playwright_regulon.py: a marker-size check asserting sizes fall in [8,26]
looked meaningful but was vacuous, because _magnitudeMarkerSizes turns every
non-finite value into 0 and then returns all-minimum 8px sizes - so a wrong
column key sailed through it. The checks that replaced it compare the
plotted values against the table they claim to come from, and were confirmed
to go red when the column key was deliberately broken.

MULTI-OMICS MOFA / FACTOR-BASED INTEGRATION (info/README-MultiOmics-MOFA.txt).
Run with the venv of requirements.txt (C:/Users/cliente/omiscera-venv): mofapy2 needs a working numpy / scipy.
Run with the server on port 8009 only for the live and browser files.

test_multiomics_mofa.py         50 checks. The engine (Multiomics_mofa), offline, on planted 3-omics data (48
                                samples; a shared factor, a transcript + protein factor, a metabolite-only
                                factor; missing protein cells). Preparation rules; refusals (< 6 samples, bad
                                settings, nothing usable); the three planted factors recovered (|r| 0.99); the
                                R2 pattern (shared vs specific); order by variance and a fixed sign; identical
                                results for the same seed and a different seed finding the same factors;
                                variance explained equal to an independent re-implementation of mofapy2's
                                definition (1e-14); stability (real factors ~1.0, injected noise factors
                                flagged); warnings; group tests against scipy, by-hand Kruskal-Wallis and the
                                BH of Genomics_module; the pathway test of the weights against scipy's
                                asymptotic Mann-Whitney for all 62 pathways with ties; planted pathways found in
                                the right direction; absolute mode; features without an identifier left out;
                                NULL CALIBRATION (random weights x random pathways: KS 0.02, 5% below 0.05);
                                collapse of several features of one gene.
test_multiomics_mofa_session.py 77 checks. The routes and jobs in API.py, offline with a stubbed resolver and
                                stubbed pathway builders: empty state; nothing allowed before the tables are
                                matched; every setting range; group sources (a file column, a sample sheet,
                                bad sheets); the run on the REAL example files (planted factors recovered, strict
                                JSON, 25 top weights are the 25 largest by |weight| of both signs, group tests);
                                groups change nothing in the fit; downloads; the enrichment job (collapse by gene,
                                the 'unavailable' path, the real builders called with the right inputs); the
                                fingerprint (same data keeps the result, a changed value or a renamed sample
                                clears it); Joint's selection keeps MOFA's state; the 6-sample floor.
                                Also: test_multiomics_correlation.py (110) gained the per-sample label columns
                                and test_multiomics_retriever.py (24) the Entrez ids / accessions kept for MOFA.
live_e2e_multiomics_mofa.py     34 checks. Real HTTP end-to-end (port 8009, real resolvers, real mofapy2, real
                                Reactome). From the raw CSV files with no Omiscera code: the three planted
                                factors recovered; the total variance explained of every omics recomputed from
                                the factor values, the downloaded weights and the raw tables (agrees to 0.5
                                point); every factor's group p-value against scipy; only the treatment factor
                                separates the groups; the top weights equal the downloaded table; the
                                pathway enrichment puts Glycolysis first; identical factors on a second run.
playwright_multiomics_mofa.py   55 checks. Real browser (port 8009): the fifth box is live; the upload screen
                                runs in MOFA mode (no correlation setup, a Continue button); the setup offers the
                                metabolite file's Group column; the three planted factors recovered; every
                                figure has the tools row, a Size menu (150% is 1.5 x wider) and Download figure;
                                the heatmap cells equal the table; the scatter has one trace per group; only the
                                treatment factor differs between groups; the weights bars are the 25 largest and
                                each has its own row (two features of one gene used to stack on one row - found
                                by looking at the screenshot); a real pathway enrichment; downloads; choosing the
                                box again keeps the matched tables; uploading again clears the match and the MOFA
                                result; the absolute-weights test is really the one run; no page error. Set MOFA_SHOT_DIR to save a screenshot of each pane.
generate_mofa_examples.py       Regenerates Data3/example_mofa_*.csv and tests/fixtures/mofa_example_truth.csv
                                (synthetic, seeded, planted structure; needs the network once for gene symbols).
                                Not a test.

MULTI-OMICS DIABLO / MULTIBLOCK INTEGRATION (info/README-MultiOmics-DIABLO.txt).
Run with the venv of requirements.txt (C:/Users/cliente/omiscera-venv). Run with the server on port 8009 only
for the live and browser files. Regenerating tests/fixtures/diablo_reference/out/*.csv needs a local R + mixOmics
(tests/r_reference/README.txt) - the offline test itself needs neither.

test_multiomics_diablo.py       52 checks. The engine (Multiomics_diablo), offline. Agreement with a REAL mixOmics
                                run on a small fixed dataset (tests/fixtures/diablo_reference/, produced once by
                                tests/r_reference/diablo_reference.R): training variates, loadings, selected-
                                feature counts, AVE (including "undefined for a block with missing values", a
                                mixOmics quirk replicated on purpose), predicted test-set scores, and every one of
                                mixOmics' three classification distances plus majority and weighted vote - all
                                equal mixOmics to float64 precision (1e-14 to 1e-15), not a looser tolerance.
                                prepare_blocks / check_groups rules; planted recovery on the REAL 36-sample example
                                (no R needed - the planted treatment genes/metabolites are selected, not the
                                unrelated cholesterol/lipid signal); NO LEAKAGE in the cross-validation, checked
                                directly (mean CV error over 12 independent pure-noise datasets is close to the
                                0.5 chance level); the permutation test's internal consistency (p = 1/(n+1) for an
                                impossible BER of 0 or 1); soft_threshold_keep's tied-value rule (a tied pair at
                                the keep/drop boundary is kept or dropped TOGETHER, R's ties.method="max");
                                weighted_vote actually uses the per-block weights (not a disguised majority vote);
                                a tied (None) prediction always counts as an error, in its TRUE class's total (not
                                excluded from both, which would wrongly look like a free pass);
                                selected_feature_enrichment's Fisher exact p-value against scipy, with the
                                restricted-universe convention the rest of the app's ORA uses.
test_multiomics_diablo_session.py 60 checks. The routes and jobs in API.py, offline with a stubbed resolver and
                                stubbed pathway builders: empty state; nothing allowed before the tables are
                                matched; a group column REQUIRED (unlike MOFA, where groups are optional); group
                                sources (a file column, a sample sheet); every setting range, including a fold
                                count bigger than the smallest group; the run on the REAL example files (strict
                                JSON - AVE's NaN for a block with missing values must serialise as null, the bug
                                found live before this test existed; component search, keepX search, a small
                                cross-validated error, a significant permutation p-value); downloads; the
                                fingerprint (a changed table clears the config, result, enrichment, cached pathway
                                sets and sheet); Joint's selection keeps DIABLO's state; the enrichment job (built
                                once per result, the real builders called with the right inputs, the 'unavailable'
                                path).
live_e2e_multiomics_diablo.py   35 checks. Real HTTP end-to-end (port 8009, real resolvers, the real from-scratch
                                DIABLO implementation, real Reactome). From the raw CSV files with no Omiscera
                                code: the planted treatment genes/proteins/metabolites are selected and correlate
                                strongly with the groups, NOT the real-but-unrelated cholesterol/lipid signal; the
                                downloaded weights and components are mutually consistent (association check); the
                                pathway enrichment puts Glycolysis/glucose metabolism/TCA cycle first, its top
                                p-value equals scipy's Fisher exact; an identical re-run gives an identical
                                cross-validated error and keepX (the seed is fixed) and correctly KEEPS the earlier
                                enrichment (still valid for the identical selection) while reconfiguring - even to
                                the same settings - correctly drops it.
playwright_multiomics_diablo.py 34 checks. Real browser (port 8009): the sixth (last) box is live; the screen says
                                DIABLO is supervised; the setup requires a group column; every result pane
                                (Overview, Components - a scatter with one marker shape per omics and colour per
                                group, so the omics' agreement is visible, and a by-group box plot; Selected
                                features; Performance - the component-count search and the permutation null
                                histogram, with the observed-BER line and label correctly inside the plot, a
                                rendering bug found by looking at the first screenshot and fixed; Pathway
                                enrichment with a real run; Data used; Parameters) has the figure tools / Size /
                                Download figure; downloads; choosing the box again keeps the matched tables;
                                uploading again clears the match and the result; no page error. Set
                                DIABLO_SHOT_DIR to save a screenshot of each pane.
generate_diablo_r_fixtures.py   Regenerates tests/fixtures/diablo_reference/{block_*.csv, group.csv,
                                train/test_samples.csv} (the small, fixed input to the mixOmics agreement check;
                                offline, reuses the MOFA example tables). Not a test; see tests/r_reference/
                                README.txt for the R half (not shipped with the repo).
