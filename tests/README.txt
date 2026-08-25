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

test_transcriptomics_network.py 44 checks. Unit tests for
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
                                both confirmed to go red).

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

test_proteomics_ora.py          21 checks. Proteomics_module.run_ora's own
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
                                183 terms tested instead of 2 - the defect
                                Genomics_module._run_fisher_for_namespace
                                STILL has, patched only locally here); and
                                (2) the universe must be the QUANTIFIED
                                proteins, not every protein the database
                                annotates - exposed with a fixture whose
                                terms deliberately contain
                                annotated-but-never-measured accessions
                                (pre-fix: K=25 instead of the correct K=5).
                                The original fixture was blind to (2), since
                                all of its genes were already quantified.


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

live_e2e_transcriptomics_network.py  22 checks. Full real HTTP round trip
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
playwright_transcriptomics_network.py  37 checks. The FIRST real-browser
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
                                appear at all.
playwright_proteomics_ora.py    40 checks. The FIRST real-browser coverage of
                                the Proteomics module (Milestone 1: the
                                foundation plus Over-Representation
                                Analysis). Drives the real flow end to end -
                                tile -> sub-tile gating (only ORA enabled;
                                Quantitative Enrichment / PPI Network / PTM
                                correctly disabled) -> organism -> upload ->
                                REAL UniProt ID Mapping resolve of all 3500
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

These five (playwright_regulon.py, playwright_genomics_results_cards.py,
playwright_msigdb.py, playwright_transcriptomics_network.py,
playwright_proteomics_ora.py) make real network calls (mygene.info, TRRUST,
ChEA3, MSigDB, UniProt, STRING, GO/QuickGO, InterPro, Complex Portal,
DISEASES) and are not mocked, so they need internet and take a few
minutes. test_msigdb_summary_figure.py, test_gsea_null_histogram_range.py,
and test_organism_allowlist_mirror.py are the exceptions - all three need
the server only to load app.js, and make no external calls.


---------------------------------------------------------------------------
RUNNING EVERYTHING
---------------------------------------------------------------------------
There is no runner script - each file is independent by design, so a failure
in one never hides another. From the repo root, with the server up:

    Get-ChildItem tests\*.py | ForEach-Object { python $_.FullName }

Expected totals as of the last full run, in file-listing order: 117, 22, 16,
4, 15, 5, 16, 10, 24, 58, 44, 49, 52, 8, 31, 21, 21, 15, 12, 12, 22, 22, 16,
12, 49, 90, 13, 37, 40, 7, 4, 15, 68 - all passing (947 checks across 33 files).

A NOTE ON WRITING NEW CHECKS HERE. Prefer assertions that can actually fail
for the reason you care about. A concrete example lives in
playwright_regulon.py: a marker-size check asserting sizes fall in [8,26]
looked meaningful but was vacuous, because _magnitudeMarkerSizes turns every
non-finite value into 0 and then returns all-minimum 8px sizes - so a wrong
column key sailed through it. The checks that replaced it compare the
plotted values against the table they claim to come from, and were confirmed
to go red when the column key was deliberately broken.
