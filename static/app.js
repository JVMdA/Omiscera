const PUBLIC_PAGES = ["Home", "Information", "Contact", "References", "Citation"];
const DISPLAY_COLUMNS = ["Name", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula", "Mapping_status"];

const STATE = {
    view: "public",
    publicPage: "Home",
    workflowStep: 1,
    aaPage: "overview",
    genomicsModuleName: "",
    genomicsIdentifierType: "",
    transcriptomicsModuleName: "",
    mc: {
        inputTable: null,
        uploadedFilename: "",
        uploadDiagnostics: null,
        mappedTable: null,
        mappedCounts: null,
        databases: [],
        selectedRanks: 1,
        ranks: {},
        selections: [],
        analysisResults: null,
        _contractWarning: null,
    },
    ora: null,
    transcriptomics: null,
    proteomics: null,
    proteomicsModuleName: "",
};

function freshOraState() {
    return {
        identifierType: "",
        taxonId: "",
        taxonName: "",
        geneTable: null,
        uploadedFilename: "",
        duplicatesRemoved: 0,
        seqInfo: null,
        resolved: null,
        ambiguousChoices: {},
        finalTable: null,
        blast: { selection: {} },
        // ORA statistical analysis (downstream of finalTable; works the same
        // regardless of which identifierType produced finalTable).
        analysisConfig: { databases: [], perDatabase: {} },
        analysisResults: null,
        analysisRunLocked: false,
        // Gene-set Burden Test - shares every field above (upload/resolve/
        // organism selection); only the analysis stage below is separate.
        scoreColumn: "",
        scoreDirection: "higher",
        scoreSummary: null,
        burdenConfig: { databases: [], perDatabase: {} },
        burdenResults: null,
        burdenRunLocked: false,
        // Network Propagation / Module Detection - shares retrieval; organism
        // comes from the STRING species selector (stored in taxonId/taxonName).
        netpropConfig: {
            moduleDetection: true, confidence: 700, fdrAlpha: 0.05,
            minModuleSize: 3, nPermutations: 1000,
        },
        netpropResults: null,
        netpropRunLocked: false,
        // Gene Function Annotation - shares retrieval (table mode). Single-gene
        // mode also writes finalTable client-side (mirroring what the server's
        // table-upload flow produces), so renderContinueToAnalysisButton and
        // the analysis screens below work identically for both modes.
        annotationSingleIdentifier: "",
        annotationSingleResult: null,
        annotationConfig: { includeIea: true },
        annotationResults: null,
        annotationRunLocked: false,
    };
}

function freshProteomicsState() {
    return {
        taxonId: "", taxonName: "",
        deTable: null, uploadedFilename: "", columnMapping: null,
        resolved: null,
        oraConfig: null, oraResults: null, oraRunLocked: false,
        oraSigAlpha: 0.05,
        customDb: null, background: null,
    };
}

function freshTranscriptomicsState() {
    return {
        identifierType: "",
        taxonId: "",
        taxonName: "",
        // Set when the organism is changed after a resolved table already
        // existed for a different one (see the organism selector's change
        // handler); shown as a warning banner until the next resolve
        // succeeds. See _clearTranscriptomicsResolutionDependentState.
        _organismChangedNotice: null,
        deTable: null,
        uploadedFilename: "",
        columnMapping: null,
        resolved: null,
        // Over-Representation Analysis (first analysis tile). Shares retrieval
        // (resolved) above; only this stage is ORA-specific.
        oraConfig: null,
        oraResults: null,
        oraRunLocked: false,
        // Display-only significance threshold for the results screen's charts/
        // row-highlighting (mirrors Genomics' per-database sigThreshold) - never
        // sent to the server, adjustable without a re-run.
        oraSigAlpha: 0.05,
        // Summary of the uploaded custom gene-set database (filename, n_terms,
        // universe_size), for display only - the parsed terms themselves live
        // server-side in the session.
        customDb: null,
        // Summary of the uploaded background/reference gene list (filename,
        // n_genes), for display only - applies globally, unlike Genomics'
        // per-database background.
        background: null,
        // Gene Set Enrichment Analysis (second analysis tile). Shares
        // retrieval AND the customDb/background widgets above with ORA
        // (same server-side session keys) - only this stage is GSEA-specific.
        gseaConfig: null,
        gseaResults: null,
        gseaRunLocked: false,
        // Display-only significance threshold for the results screen (mirrors
        // oraSigAlpha above) - never sent to the server. 0.25 is GSEA's own
        // convention default (looser than ORA's 0.05), since GSEA's
        // permutation FDR is already a stricter correction on its own.
        gseaSigAlpha: 0.25,
        // Transcription Factor / Regulon (third analysis tile). Shares
        // retrieval above; a separate config/results pair per analysis
        // (target enrichment vs. activity inference), since a run can select
        // either or both (server: Transcriptomics_regulon.py).
        regulonConfig: null,
        regulonResults: null,
        regulonRunLocked: false,
        regulonSigAlpha: 0.05,
        // Network / Module Detection (fourth analysis tile). Also shares the
        // retrieval above; seeds the STRING interaction network with the DEGs
        // (server: Transcriptomics_network.py).
        networkConfig: null,
        networkResults: null,
        networkRunLocked: false,
        networkSigAlpha: 0.05,
    };
}

// Clears every piece of Transcriptomics state that depends on the resolved
// table - shared by the organism-change handler and the upload handler, so
// the two "something upstream changed, discard downstream analyses" call
// sites can't drift apart. They used to: the upload handler cleared
// resolved/customDb/background, while the SERVER's own upload reset
// (_fresh_transcriptomics_state() in API.py) replaces the entire session
// sub-dict at once - so after a new upload the client could still show ORA/
// GSEA/regulon results the server had already discarded underneath it.
//
// Deliberately does NOT touch identifierType, deTable/uploadedFilename,
// customDb or background - those are organism-independent and the caller's
// own business (the upload handler clears customDb/background itself,
// deliberately, as part of "new file = fresh start"; an organism change
// alone should NOT discard a still-relevant custom gene-set upload).
function _clearTranscriptomicsResolutionDependentState() {
    const tx = STATE.transcriptomics;
    tx.resolved = null;
    tx.oraConfig = null; tx.oraResults = null; tx.oraRunLocked = false;
    tx.gseaConfig = null; tx.gseaResults = null; tx.gseaRunLocked = false;
    tx.regulonConfig = null; tx.regulonResults = null; tx.regulonRunLocked = false;
    tx.networkConfig = null; tx.networkResults = null; tx.networkRunLocked = false;
    // Draft config (selected databases/collections/thresholds) is lazily
    // created per config screen (STATE.transcriptomics._oraDraft etc.) and
    // is organism-gated (e.g. MSigDB collection availability) - deleting it
    // forces a fresh, correctly-gated draft next time that screen renders,
    // rather than carrying over picks that may no longer be valid.
    delete tx._oraDraft;
    delete tx._gseaDraft;
    delete tx._regulonDraft;
    delete tx._networkDraft;
}

// Coverage per regulon database (see Transcriptomics_regulon.py's module
// docstring) - used to grey out a database's checkbox for an organism it has
// no data for, the same idiom TRANSCRIPTOMICS_ORA_DATABASES uses for MSigDB
// ("human/mouse only", plus a per-collection check - see
// _msigdbCollectionAvailableForTaxon).
// Below this many resolved genes, a TF/Regulon run will realistically return
// nothing: a TF is only tested when "Minimum regulon size" (default 5) of its
// targets are present, and regulons are a small slice of the transcriptome.
// Used only to warn up front - it never blocks a run.
const REGULON_MIN_GENES_FOR_USEFUL_RUN = 500;

const TRANSCRIPTOMICS_REGULON_DATABASES = ["TRRUST", "DoRothEA", "CollecTRI", "ChEA3"];
const TRANSCRIPTOMICS_REGULON_TAXA = { "9606": true, "10090": true };  // human + mouse
const TRANSCRIPTOMICS_REGULON_DB_COVERAGE = {
    "TRRUST": "human & mouse only", "DoRothEA": "human & mouse only",
    "CollecTRI": "human & mouse only", "ChEA3": "human only",
};
function _regulonDbAvailableForTaxon(dbName, taxonId) {
    if (dbName === "ChEA3") return taxonId === "9606";
    return !!TRANSCRIPTOMICS_REGULON_TAXA[taxonId];
}

const TRANSCRIPTOMICS_ID_LABELS = {
    gene_symbol: "Gene Symbol",
    ensembl_gene: "Ensembl Gene ID",
    ensembl_transcript: "Ensembl Transcript ID",
    entrez: "Entrez Gene ID",
};
// The name the browser SAVES each example as. MUST stay in sync with API.py's
// TRANSCRIPTOMICS_EXAMPLE_BY_IDENTIFIER_TYPE, which decides what is actually
// served - downloadFile() sets a.download from this map, so a mismatch means
// the user saves real file X under the name of unrelated file Y. That is
// exactly what happened here: gene_symbol served the 3000-gene de_full file
// but saved it as "example_transcriptomics_symbol.csv", the name of a real,
// different 20-gene file. A consistency test now parses both files and
// asserts these two maps are identical.
const TRANSCRIPTOMICS_EXAMPLE_FILENAMES = {
    gene_symbol: "example_transcriptomics_de_full.csv",
    ensembl_gene: "example_transcriptomics_ensembl_full.csv",
    ensembl_transcript: "example_transcriptomics_transcript.csv",
    entrez: "example_transcriptomics_entrez_full.csv",
};

// Databases enabled in the ORA config screen. Extended as each database's
// universe builder is implemented and verified (Reactome, then GWAS Catalog).
const ORA_CUSTOM_DB_LABEL = "Custom database";
const ORA_AVAILABLE_DATABASES = ["GO", "Reactome", "GWAS Catalog", ORA_CUSTOM_DB_LABEL];
const ORA_DATABASE_CAPTIONS = {
    "GO": "Gene Ontology. Produces three independent result tables: Biological Process, Molecular Function, Cellular Component. Each is tested and FDR-corrected separately.",
    "Reactome": "Pathway over-representation using UniProt's Reactome cross-references, tested with our own Fisher's exact + BH-FDR pipeline.",
    "GWAS Catalog": "Trait associations from the EBI GWAS Catalog. Works only with Gene Symbol input - matches your original gene symbols directly against the catalog. Coverage is overwhelmingly human; non-human species selections will likely return few or zero results.",
    "Custom database": "Upload your own gene sets (GMT or a tidy term/gene table) and test against them with the same Fisher's exact + BH-FDR pipeline.",
};
const ORA_GO_NAMESPACE_LABELS = {
    "biological_process": "Biological Process",
    "molecular_function": "Molecular Function",
    "cellular_component": "Cellular Component",
};

// Gene-set Burden Test: same three families as ORA minus GWAS Catalog (an
// unconventional/circular unit for a variant-burden test) and no background
// upload option (not applicable to a competitive rank test).
const BURDEN_AVAILABLE_DATABASES = ["GO", "Reactome", ORA_CUSTOM_DB_LABEL];
const BURDEN_DATABASE_CAPTIONS = {
    "GO": "Gene Ontology. Produces three independent result tables: Biological Process, Molecular Function, Cellular Component. Each is tested and FDR-corrected separately.",
    "Reactome": "Pathway burden test using UniProt's Reactome cross-references, tested with a competitive Wilcoxon rank-sum + BH-FDR pipeline.",
    "Custom database": "Upload your own gene sets (GMT or a tidy term/gene table) and test against them with the same competitive Wilcoxon + BH-FDR pipeline.",
};
const ORA_DOT_COLOR_SCALES = [
    ["YlOrRd", "YlOrRd", [[0, "#ffffcc"], [0.35, "#fed976"], [0.7, "#fd8d3c"], [1, "#800026"]], false],
    ["YlGnBu", "YlGnBu", [[0, "#ffffd9"], [0.35, "#c7e9b4"], [0.7, "#41b6c4"], [1, "#081d58"]], false],
    ["Portland", "Portland", "Portland", false],
    ["Purples", "Purples", "Purples", false],
    ["Greens", "Greens", [[0, "#f7fcf5"], [0.35, "#c7e9c0"], [0.7, "#41ab5d"], [1, "#00441b"]], false],
];
const ORA_BAR_COLORS = [
    ["Blue", "#1B639C"],
    ["Teal", "#158A8A"],
    ["Green", "#2E8B57"],
    ["Purple", "#7057A3"],
    ["Red", "#B84A52"],
    ["Gray", "#6C7A89"],
];

function resetMcState() {
    STATE.mc = {
        inputTable: null,
        uploadedFilename: "",
        uploadDiagnostics: null,
        mappedTable: null,
        mappedCounts: null,
        databases: STATE.mc.databases,
        selectedRanks: 1,
        ranks: {},
        selections: [],
        analysisResults: null,
        _contractWarning: STATE.mc._contractWarning,
    };
    STATE.workflowStep = 1;
}

function getSessionId() {
    let id = localStorage.getItem("omiscera_session_id");
    if (!id) {
        id = (crypto.randomUUID ? crypto.randomUUID() : String(Date.now()) + Math.random().toString(16).slice(2));
        localStorage.setItem("omiscera_session_id", id);
    }
    return id;
}

async function apiFetch(path, options) {
    options = options || {};
    const headers = Object.assign({ "X-Session-Id": getSessionId() }, options.headers || {});
    const resp = await fetch(path, Object.assign({}, options, { headers }));
    if (!resp.ok) {
        let detail = resp.statusText;
        try {
            const data = await resp.clone().json();
            if (data && data.detail) detail = data.detail;
        } catch (e) { /* response was not JSON */ }
        throw new Error(detail);
    }
    return resp;
}

async function apiJson(path, options) {
    const resp = await apiFetch(path, options);
    return resp.json();
}

async function apiPostJson(path, body) {
    return apiJson(path, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
    });
}

async function apiUpload(path, formData) {
    return apiJson(path, { method: "POST", body: formData });
}

async function downloadFile(path, filename) {
    try {
        const resp = await apiFetch(path);
        const blob = await resp.blob();
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
    } catch (err) {
        alert("Download failed: " + err.message);
    }
}

function pollJob(jobId, callbacks) {
    const onProgress = callbacks.onProgress || function () {};
    const onDone = callbacks.onDone || function () {};
    const onError = callbacks.onError || function () {};
    const interval = setInterval(async () => {
        try {
            const job = await apiJson("/api/jobs/" + jobId);
            onProgress(job);
            if (job.status === "done") {
                clearInterval(interval);
                onDone(job.result);
            } else if (job.status === "error") {
                clearInterval(interval);
                onError(job.error || "Job failed.");
            }
        } catch (err) {
            clearInterval(interval);
            onError(err.message);
        }
    }, 1100);
}

function el(tag, props, children) {
    props = props || {};
    const node = document.createElement(tag);
    Object.keys(props).forEach((key) => {
        const value = props[key];
        if (value === null || value === undefined) return;
        if (key === "class") node.className = value;
        else if (key === "style" && typeof value === "object") Object.assign(node.style, value);
        else if (key.indexOf("on") === 0 && typeof value === "function") node.addEventListener(key.slice(2), value);
        else node.setAttribute(key, value);
    });
    const childList = children === undefined ? [] : (Array.isArray(children) ? children : [children]);
    childList.forEach((child) => {
        if (child === null || child === undefined) return;
        if (typeof child === "string" || typeof child === "number") node.appendChild(document.createTextNode(String(child)));
        else node.appendChild(child);
    });
    return node;
}

function alertBox(kind, text) {
    return el("div", { class: "alert alert-" + kind }, text);
}

function truncate(value, maxLen) {
    const text = value === null || value === undefined ? "" : String(value);
    return text.length > maxLen ? text.slice(0, maxLen) + "..." : text;
}

function fileInputWidget(opts) {
    opts = opts || {};
    const input = el("input", { type: "file", class: "visually-hidden-file-input" });
    if (opts.accept) input.setAttribute("accept", opts.accept);
    const filenameSpan = el("span", { class: "file-name-display" }, opts.initialFilename || "No file selected");
    const button = el("button", { type: "button", class: "file-choose-btn" }, opts.buttonLabel || "Choose file");
    button.addEventListener("click", () => input.click());
    input.addEventListener("change", () => {
        const file = input.files && input.files[0] ? input.files[0] : null;
        filenameSpan.textContent = file ? file.name : "No file selected";
        if (opts.onChange) opts.onChange(file);
    });
    return el("div", { class: "file-drop" }, [button, filenameSpan, input]);
}

// Columns across the app that hold a comma-joined gene list, sometimes
// hundreds of symbols long (MSigDB Hallmark sets routinely overlap 100+
// DEGs). Without special handling, one such cell wrapping inside a narrow
// column inflates its WHOLE ROW to thousands of pixels tall - and since
// table cells default to vertical-align:middle, every other (short) cell in
// that row gets centered deep inside that height, off the bottom of the
// results table's scroll viewport. The table then looks completely empty:
// the data is there, just rendered ~1500px below what's visible. Both
// halves of the fix below address this: capping these specific cells stops
// the row-height explosion at the source, and vertical-align:top is kept as
// a second, independent safeguard for any other column that turns out to
// hold unexpectedly long content.
// "Overlapping Proteins" is Proteomics' own name for the same kind of long
// comma-joined member list (Proteomics_module.OVERLAP_COLUMN) - it MUST be
// listed here too, or that module's result tables reintroduce the exact
// row-height explosion this clamp exists to prevent (see
// tests/test_table_long_list_rendering.py).
const TABLE_LONG_LIST_COLUMNS = new Set([
    "Overlapping Genes", "Overlapping Proteins", "Leading Edge Genes", "Genes",
]);

// scipy p-values/FDRs (hypergeom.sf, a normal-tail p, BH correction on either)
// can genuinely underflow to a literal 0.0 for an astronomically significant
// result - too small for float64 to represent (below ~1e-308). Naively
// formatting that as toExponential(N) renders "0.000e+0", which reads exactly
// like a computation bug rather than what it is. Shared by every p-value/FDR
// column in the app (ORA, GSEA, Regulon enrichment + activity, Network hits +
// modules) so the fix can't drift out of sync between tables. buildTable
// already turns a null/undefined cell into "" before calling cellFormatFn, so
// that (not null) is the guard needed here - without it Number("") is 0,
// which is finite, and every blank cell (e.g. a recruited gene's padj) would
// wrongly render as "< 1e-300" instead of staying blank.
function _fmtSci(val, digits) {
    if (val === "" || val === null || val === undefined) return "";
    const num = Number(val);
    if (!isFinite(num)) return val;
    if (num === 0) return "< 1e-300";
    return num.toExponential(digits === undefined ? 3 : digits);
}
function _fmtFixed(val, digits) {
    if (val === "" || val === null || val === undefined) return "";
    const num = Number(val);
    return isFinite(num) ? num.toFixed(digits) : val;
}

function buildTable(payload, opts) {
    opts = opts || {};
    const wrap = el("div", { class: "table-wrap" });
    const scroller = el("div", { class: "table-scroll" });
    if (opts.maxHeight) scroller.style.maxHeight = opts.maxHeight + "px";
    const table = el("table", { class: "data-table" });
    const columns = payload && payload.columns ? payload.columns : [];
    const rows = payload && payload.rows ? payload.rows : [];

    const thead = el("thead");
    const headRow = el("tr");
    columns.forEach((col) => headRow.appendChild(el("th", {}, col)));
    thead.appendChild(headRow);
    table.appendChild(thead);

    const tbody = el("tbody");
    rows.forEach((row, idx) => {
        const cls = opts.rowClassFn ? opts.rowClassFn(row, idx) : null;
        const tr = el("tr", cls ? { class: cls } : {});
        columns.forEach((col) => {
            let val = row[col];
            if (val === null || val === undefined) val = "";
            if (opts.cellFormatFn) val = opts.cellFormatFn(col, val, row);
            const text = String(val);
            if (TABLE_LONG_LIST_COLUMNS.has(col) && text.length > 100) {
                // Clamped to 3 lines with an ellipsis; the FULL list is still
                // available via the title tooltip on hover and, uncapped, in
                // the CSV download - this is a display cap only, never a
                // data cap.
                const td = el("td", { class: "td-long-list", title: text });
                td.appendChild(el("div", { class: "td-long-list-content" }, text));
                tr.appendChild(td);
            } else {
                tr.appendChild(el("td", {}, text));
            }
        });
        tbody.appendChild(tr);
    });
    table.appendChild(tbody);
    scroller.appendChild(table);
    wrap.appendChild(scroller);
    return wrap;
}

function metricCards(items, columns) {
    const grid = el("div", { class: "metric-cards", style: { gridTemplateColumns: "repeat(" + columns + ", 1fr)" } });
    items.forEach(([label, value]) => {
        grid.appendChild(el("div", { class: "metric-card" }, [
            el("div", { class: "metric-value" }, String(value)),
            el("div", { class: "metric-label" }, label),
        ]));
    });
    return grid;
}

// Section navigation rail for results screens: a permanent right-hand map of
// every section a run actually produced, with exactly one section visible in
// the body at a time - replaces the old single flat scroll (dropdown, table
// toggle, N chart-group toggles all stacked one after another).
//
// hostEl gets EMPTIED and given the rail layout; callers then call
// nav.pane(title, opts) for each section and append their content into the
// returned <section>. opts.group nests the pane under a non-clickable
// sub-heading in the rail (used for GSEA's six chart groups under "Charts").
//
// opts.stateKey, if given, is a key into STATE.transcriptomics/ora/etc under
// which the active pane's title is remembered across a renderPanel() rebuild
// (changing the significance threshold or the database select tears down and
// rebuilds the whole panel) - without this, nudging a threshold would always
// throw the user back to the first section.
function makeResultNav(hostEl, opts) {
    opts = opts || {};
    hostEl.innerHTML = "";
    hostEl.classList.add("result-nav-layout");

    const body = el("div", { class: "result-nav-body" });
    const rail = el("nav", { class: "result-nav" });
    hostEl.appendChild(body);
    hostEl.appendChild(rail);

    const panes = [];   // { title, section, btn }
    let activeTitle = null;
    let currentGroupWrap = null;
    let currentGroupTitle = null;

    function rememberActive() {
        if (opts.stateKey) STATE._resultNavActive = STATE._resultNavActive || {};
        if (opts.stateKey) STATE._resultNavActive[opts.stateKey] = activeTitle;
    }

    function activate(title) {
        const found = panes.find(p => p.title === title);
        if (!found) return;
        activeTitle = title;
        panes.forEach(p => {
            const isActive = p === found;
            p.section.style.display = isActive ? "" : "none";
            p.btn.classList.toggle("active", isActive);
        });
        rememberActive();
        // A Plotly chart laid out while its pane is display:none renders at
        // zero width and stays that way once shown again. Plotly.Plots.resize()
        // alone is NOT enough here - confirmed live it rescales the plot's
        // overall size but leaves stale, wrongly-computed things like legend
        // height/position exactly as they were miscalculated at the original
        // (zero-width) draw, e.g. a horizontal legend pinned ~20px too high,
        // overlapping the xaxis title. _newOraPlot stashes a full redraw
        // closure on the div for exactly this: prefer it (a genuinely fresh
        // Plotly.newPlot at the real, now-correct container size) and only
        // fall back to resize() for any chart that predates/bypasses it.
        if (typeof Plotly !== "undefined") {
            found.section.querySelectorAll(".chart-container").forEach((div) => {
                if (!div.data) return;
                try {
                    if (typeof div._omiscerRedraw === "function") {
                        div._omiscerRedraw();
                    } else {
                        Plotly.Plots.resize(div);
                    }
                } catch (e) { /* not a plot */ }
            });
        }
    }

    function pane(title, paneOpts) {
        paneOpts = paneOpts || {};
        const section = el("section", { class: "result-pane" });
        // The section title lives in the rail, but the content area still needs
        // one - without it a pane opens on a caption, some buttons and a chart
        // with no heading anywhere. Header text is the SAME string as the rail
        // entry, so the two can never drift.
        const header = el("header", { class: "result-pane-header" },
            el("h2", { class: "result-pane-title" }, title));
        section.appendChild(header);

        // Callers get an inner body div, NOT the section. Several call sites do
        // `pane.insertBefore(note, pane.firstChild)` to put a threshold note
        // above their controls; if they were handed the section, firstChild
        // would be the header and those notes would land above the title.
        // Returning the body keeps firstChild meaning "first thing the caller
        // appended", so every existing call site works unchanged.
        const paneBody = el("div", { class: "result-pane-body" });
        section.appendChild(paneBody);
        body.appendChild(section);

        if (paneOpts.group !== currentGroupTitle) {
            currentGroupTitle = paneOpts.group;
            currentGroupWrap = null;
            if (paneOpts.group) {
                rail.appendChild(el("div", { class: "result-nav-group-label" }, paneOpts.group));
                currentGroupWrap = el("div", { class: "result-nav-group" });
                rail.appendChild(currentGroupWrap);
            }
        }
        const btn = el("button", { class: "result-nav-btn" + (paneOpts.group ? " result-nav-btn-nested" : ""), type: "button" }, title);
        btn.addEventListener("click", () => activate(title));
        (currentGroupWrap || rail).appendChild(btn);

        panes.push({ title, section, btn });
        return paneBody;
    }

    return {
        pane,
        activate,
        current: () => activeTitle,
        finish: () => {
            if (panes.length === 0) return;
            const remembered = opts.stateKey && STATE._resultNavActive && STATE._resultNavActive[opts.stateKey];
            activate(remembered && panes.some(p => p.title === remembered) ? remembered : panes[0].title);
        },
    };
}

// Chart queue: populated during render, flushed via rAF after DOM update.
let _chartQueue = [];

function scheduleChart(id, buildFn) {
    _chartQueue.push({ id, buildFn });
}

function render() {
    _chartQueue = [];
    const app = document.getElementById("app");
    app.innerHTML = "";
    app.appendChild(renderTopNav());
    if (STATE.view === "metabolic_connection") app.appendChild(renderMetabolicConnection());
    else if (STATE.view === "annotation_analytics") app.appendChild(renderAnnotationAnalytics());
    else if (STATE.publicPage === "Home") app.appendChild(renderHome());
    else app.appendChild(renderPlaceholder(STATE.publicPage));
    if (_chartQueue.length > 0) {
        const queue = _chartQueue.slice();
        requestAnimationFrame(() => {
            queue.forEach(({ id, buildFn }) => {
                const div = document.getElementById(id);
                if (div) {
                    try { buildFn(div); } catch (e) { /* Plotly not ready or data missing */ }
                }
            });
        });
    }
}

function setPage(pageName) {
    if (PUBLIC_PAGES.indexOf(pageName) !== -1) {
        if (STATE.view === "metabolic_connection") resetMcState();
        STATE.publicPage = pageName;
        STATE.view = "public";
    } else if (pageName === "Annotation Analytics") {
        if (STATE.view === "metabolic_connection") resetMcState();
        STATE.view = "annotation_analytics";
        STATE.aaPage = "overview";
    }
    render();
}

function startMetabolicConnectionWorkflow() {
    resetMcState();
    STATE.view = "metabolic_connection";
    STATE.workflowStep = 1;
    render();
}

function setWorkflowStep(step) {
    STATE.workflowStep = step;
    STATE.view = "metabolic_connection";
    render();
}

function renderTopNav() {
    const bar = el("div", { class: "top-nav-bar" });
    PUBLIC_PAGES.forEach((page) => {
        const isActive = STATE.view === "public" && STATE.publicPage === page;
        bar.appendChild(el("button", {
            class: "nav-btn" + (isActive ? " active" : ""),
            onclick: () => setPage(page),
        }, page));
    });
    return bar;
}

function renderHome() {
    const wrap = el("div", { class: "home-screen" });
    wrap.appendChild(el("img", { class: "home-logo", src: "/api/assets/logo.png", alt: "Omiscera" }));
    wrap.appendChild(el("div", { class: "home-heading" }, [
        el("div", { class: "home-subtitle" }, "A Platform for Omics Evidence, Visualization and Analytics."),
    ]));
    const buttons = el("div", { class: "home-buttons" });
    buttons.appendChild(el("button", { class: "big-btn", onclick: startMetabolicConnectionWorkflow }, "Metabolic Connection"));
    buttons.appendChild(el("button", { class: "big-btn", onclick: () => setPage("Annotation Analytics") }, "Annotation Analytics"));
    wrap.appendChild(buttons);
    return wrap;
}

function renderPlaceholder(title) {
    const wrap = el("div", {});
    wrap.appendChild(el("h2", {}, title));
    wrap.appendChild(el("p", {}, "This section is not implemented in the first desktop prototype. Please use the Metabolic Connection workflow."));
    return wrap;
}

function renderMetabolicConnection() {
    const wrap = el("div", { class: "mc-module-screen" });
    wrap.appendChild(el("div", { class: "aa-header" }, [
        el("img", { src: "/api/assets/logo_MC.png", alt: "Metabolic Connection" }),
        el("div", { class: "title" }, "Metabolic Connection"),
    ]));
    if (STATE.workflowStep === 1) wrap.appendChild(renderMcStep1());
    else if (STATE.workflowStep === 2) wrap.appendChild(renderMcStep2());
    else if (STATE.workflowStep === 3) wrap.appendChild(renderMcStep3());
    else if (STATE.workflowStep === 4) wrap.appendChild(renderMcStep4());
    return wrap;
}

function renderMcStep1() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "step-subheader" }, "Step 1 - Upload compound table"));

    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.mc.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            try {
                const data = await apiUpload("/api/mc/upload", formData);
                STATE.mc.inputTable = data.table;
                STATE.mc.uploadedFilename = file.name;
                STATE.mc.uploadDiagnostics = data;
                STATE.mc.mappedTable = null;
                STATE.mc.mappedCounts = null;
                STATE.mc.analysisResults = null;
                render();
            } catch (err) {
                alert("Unable to read CSV file: " + err.message);
            }
        },
    }));
    wrap.appendChild(el("div", { class: "caption" }, "Required column: Name. Optional columns: PubChem ID, SMILES, InChIKey, HMDB ID, Formula."));

    if (STATE.mc.inputTable) {
        const maxHeight = Math.min(900, Math.max(240, (STATE.mc.inputTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.mc.inputTable, { maxHeight }));

        const diag = STATE.mc.uploadDiagnostics;
        const details = el("details", { class: "expander" });
        details.appendChild(el("summary", {}, "Upload/import diagnostics"));
        const body = el("div", { class: "expander-body" });
        if (diag.empty_name_count > 0) {
            body.appendChild(alertBox("warning", diag.empty_name_count + " empty Name value(s) were found."));
        }
        const cards = el("div", { class: "diagnostic-cards", style: { gridTemplateColumns: "repeat(4, 1fr)" } });
        [
            ["Duplicated compound names", diag.duplicate_name_count],
            ["Duplicated rows removed", diag.duplicate_rows_removed],
            ["Rows before duplicate removal", diag.rows_before_duplicate_removal],
            ["Rows after duplicate removal", diag.rows_after_duplicate_removal],
        ].forEach(([label, value]) => {
            cards.appendChild(el("div", { class: "diagnostic-card" }, [
                el("div", { class: "diagnostic-value" }, String(value)),
                el("div", { class: "diagnostic-label" }, label),
            ]));
        });
        body.appendChild(cards);
        if (diag.duplicate_report && diag.duplicate_report.rows.length > 0) {
            body.appendChild(alertBox("warning", "Duplicated compound names were detected and removed before mapping."));
            body.appendChild(buildTable(diag.duplicate_report, {}));
        } else {
            body.appendChild(alertBox("info", "No duplicated compound names were detected."));
        }
        details.appendChild(body);
        wrap.appendChild(details);
    }

    const actions = el("div", { class: "right-actions" });
    const btnRow = el("div", { class: "btn-row" });
    const continueBtn = el("button", { class: "btn" }, "Continue to compound mapping");
    continueBtn.disabled = !STATE.mc.inputTable;
    continueBtn.addEventListener("click", () => setWorkflowStep(2));
    btnRow.appendChild(continueBtn);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);
    return wrap;
}

function renderMcStep2() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "step-subheader" }, "Step 2 - Compound mapping"));

    if (STATE.mc.inputTable && !STATE.mc.mappedTable && STATE.mc.inputTable.rows.length > 50) {
        wrap.appendChild(alertBox("warning", "Large tables may take several minutes on the first run if compounds are not cached"));
    }

    const row = el("div", { class: "row" });
    const runBtn = el("button", { class: "btn" }, "Run compound mapping");
    row.appendChild(el("div", { class: "col", style: { flex: "0 0 240px" } }, runBtn));
    const statusCol = el("div", { class: "col" }, el("div", { class: "caption" },
        STATE.mc.mappedTable ? "Compound mapping completed." : "No mapping has been run yet."));
    row.appendChild(statusCol);
    wrap.appendChild(row);

    const tableContainer = el("div", {});
    const metricsContainer = el("div", {});
    if (STATE.mc.mappedTable) {
        tableContainer.appendChild(buildTable(STATE.mc.mappedTable, {}));
        const counts = STATE.mc.mappedCounts;
        metricsContainer.appendChild(metricCards([
            ["Mapped", counts.mapped], ["Partial", counts.partial], ["Not found", counts.not_found],
        ], 3));
    } else {
        tableContainer.appendChild(buildTable({ columns: DISPLAY_COLUMNS, rows: [] }, {}));
        metricsContainer.appendChild(metricCards([["Mapped", 0], ["Partial", 0], ["Not found", 0]], 3));
    }
    wrap.appendChild(tableContainer);
    wrap.appendChild(metricsContainer);

    if (STATE.mc.mappedTable) {
        const dlBtn = el("button", { class: "btn" }, "Download mapped compounds");
        dlBtn.addEventListener("click", () => downloadFile("/api/mc/download/mapped.csv", "mapped_compounds.csv"));
        wrap.appendChild(dlBtn);
    }

    runBtn.addEventListener("click", async () => {
        if (!STATE.mc.inputTable) { alert("Please upload a CSV file before running mapping."); return; }
        runBtn.disabled = true;
        statusCol.innerHTML = "";
        statusCol.appendChild(el("div", { class: "caption" }, "Running compound mapping..."));
        try {
            const { job_id } = await apiPostJson("/api/mc/map", {});
            pollJob(job_id, {
                onDone: (result) => {
                    STATE.mc.mappedTable = result.table;
                    STATE.mc.mappedCounts = result.counts;
                    STATE.mc.analysisResults = null;
                    render();
                },
                onError: (err) => { alert("Compound mapping failed: " + err); render(); },
            });
        } catch (err) {
            alert("Compound mapping failed: " + err.message);
            render();
        }
    });

    const actions = el("div", { class: "right-actions" });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back to Upload");
    backBtn.addEventListener("click", () => { resetMcState(); setWorkflowStep(1); });
    const continueBtn = el("button", { class: "btn" }, "Continue to Rank Selection");
    continueBtn.disabled = !STATE.mc.mappedTable;
    continueBtn.addEventListener("click", () => setWorkflowStep(3));
    btnRow.appendChild(backBtn);
    btnRow.appendChild(continueBtn);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);
    return wrap;
}

function getRankState(letter) {
    if (!STATE.mc.ranks[letter]) STATE.mc.ranks[letter] = { selectedDbs: [], contexts: {}, hmdb: null, hmdbPreview: null };
    return STATE.mc.ranks[letter];
}

function buildSelectionsList() {
    const selections = [];
    Object.keys(STATE.mc.ranks).sort().forEach((letter) => {
        if (letter.charCodeAt(0) - 65 >= STATE.mc.selectedRanks) return;
        const rankState = STATE.mc.ranks[letter];
        rankState.selectedDbs.forEach((dbName) => {
            if (dbName === "HMDB Import") {
                if (rankState.hmdb) selections.push(rankState.hmdb);
            } else if (rankState.contexts[dbName]) {
                selections.push(rankState.contexts[dbName]);
            }
        });
    });
    return selections;
}

function renderMcStep3() {
    const wrap = el("div", { class: "mc-step3-screen" });
    wrap.appendChild(el("div", { class: "step-subheader" }, "Step 3 - Database and biological context selection"));
    const panel = el("div", { class: "mc-step3-panel" });
    wrap.appendChild(panel);
    ensureMcStep3Loaded(panel);
    return wrap;
}

async function ensureMcStep3Loaded(panel) {
    if (!STATE.mc.databases || STATE.mc.databases.length === 0) {
        try {
            const data = await apiJson("/api/databases");
            STATE.mc.databases = data.databases;
        } catch (err) {
            STATE.mc.databases = ["LOTUS", "Rhea + UniProt", "NPAtlas", "HMDB Import"];
        }
    }
    if (STATE.mc._contractWarning === null) {
        try {
            const contract = await apiJson("/api/catalog/validate");
            STATE.mc._contractWarning = contract.ok ? "" :
                "database_catalog.py contract validation reported issues: " + (contract.errors || []).join("; ");
        } catch (err) {
            STATE.mc._contractWarning = "database_catalog.py contract validation could not run: " + err.message;
        }
    }
    renderStep3Panel(panel);
}

function renderStep3Panel(panel) {
    panel.innerHTML = "";
    if (STATE.mc._contractWarning) panel.appendChild(alertBox("warning", STATE.mc._contractWarning));

    const rankRow = el("div", { class: "mc-rank-count-row" });
    rankRow.appendChild(el("label", { class: "field-label" }, "Number of ranks:"));
    const rankSelect = el("select", {});
    for (let i = 1; i <= 8; i++) rankSelect.appendChild(el("option", { value: i }, String(i)));
    rankSelect.value = String(STATE.mc.selectedRanks);
    rankSelect.addEventListener("change", () => {
        STATE.mc.selectedRanks = parseInt(rankSelect.value, 10);
        renderStep3Panel(panel);
    });
    rankRow.appendChild(rankSelect);
    panel.appendChild(rankRow);

    for (let i = 0; i < STATE.mc.selectedRanks; i++) {
        const rankLetter = String.fromCharCode(65 + i);
        panel.appendChild(renderRankCard(rankLetter, panel));
    }

    panel.appendChild(renderSelectionsSummary());

    const actions = el("div", { class: "right-actions" });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back to mapping review");
    backBtn.addEventListener("click", () => setWorkflowStep(2));
    const runBtn = el("button", { class: "btn" }, "Run connection analysis");
    runBtn.addEventListener("click", async () => {
        const selections = buildSelectionsList();
        if (selections.length === 0) {
            alert("Select at least one database/context before running the analysis.");
            return;
        }
        runBtn.disabled = true;
        try {
            const { job_id } = await apiPostJson("/api/mc/analyze", { selections });
            pollJob(job_id, {
                onDone: (result) => { STATE.mc.analysisResults = result; setWorkflowStep(4); },
                onError: (err) => { alert("Connection analysis failed: " + err); runBtn.disabled = false; },
            });
        } catch (err) {
            alert("Connection analysis failed: " + err.message);
            runBtn.disabled = false;
        }
    });
    btnRow.appendChild(backBtn);
    btnRow.appendChild(runBtn);
    actions.appendChild(btnRow);
    panel.appendChild(actions);
}

const DB_CONTEXT_CAPTIONS = {
    "LOTUS": "Context selection uses Wikidata taxon search as a taxonomic helper. Enter an organism name or - preferably - its Wikidata Q-code for an exact match. Note: several entities can share the same name, so a name search may surface the wrong one; when unsure, verify the taxon on Wikidata and LOTUS and paste its Q-code.",
    "Rhea + UniProt": "Context is a UniProt taxon. Enter a species name or numeric taxon ID for an exact match. Using the numeric taxon ID is recommended to avoid ambiguity - you can look it up at uniprot.org/taxonomy. Enzyme presence is checked against all UniProt entries (reviewed + unreviewed) for the chosen taxon, using exact Rhea reaction IDs.",
    "NPAtlas": "NPAtlas covers microbial natural products (bacteria + fungi). Search by genus or family name to match compounds produced by organisms within that group. You can also type a full species name for an exact species match - select the row labelled with your typed text for direct name matching.",
    "HMDB Import": "Upload an HMDB export file (no extension, CSV, or XLSX). One file per rank. Compound lookup will run in the next phase.",
};

function renderRankCard(rankLetter, panel) {
    const rankState = getRankState(rankLetter);
    const card = el("div", { class: "rank-card mc-rank-card mc-rank-card-" + rankLetter.toLowerCase() });
    card.appendChild(el("div", { class: "section-label mc-rank-title" }, "Rank " + rankLetter));

    const checkRow = el("div", { class: "db-checkbox-row mc-rank-db-row" });
    STATE.mc.databases.forEach((dbName) => {
        const label = el("label", { class: "mc-db-choice" });
        const cb = el("input", { type: "checkbox" });
        cb.checked = rankState.selectedDbs.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (rankState.selectedDbs.indexOf(dbName) === -1) rankState.selectedDbs.push(dbName);
            } else {
                rankState.selectedDbs = rankState.selectedDbs.filter((d) => d !== dbName);
                delete rankState.contexts[dbName];
                if (dbName === "HMDB Import") { rankState.hmdb = null; rankState.hmdbPreview = null; }
            }
            renderStep3Panel(panel);
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName));
        checkRow.appendChild(label);
    });
    card.appendChild(checkRow);

    rankState.selectedDbs.forEach((dbName) => {
        const block = el("div", { style: { margin: "0.6rem 0 1rem 0" } });
        block.appendChild(el("div", { class: "section-label", style: { margin: "0.4rem 0" } }, dbName));
        if (DB_CONTEXT_CAPTIONS[dbName]) block.appendChild(el("div", { class: "caption" }, DB_CONTEXT_CAPTIONS[dbName]));
        if (dbName === "HMDB Import") block.appendChild(renderHmdbImportBlock(rankLetter, rankState, panel));
        else block.appendChild(renderContextSearchWidget(rankLetter, dbName, rankState, panel));
        card.appendChild(block);
    });

    return card;
}

function renderContextSearchWidget(rankLetter, dbName, rankState, panel) {
    const wrap = el("div", {});
    const input = el("input", { class: "text-input", type: "text", placeholder: dbName + " context (organism name or taxon ID)" });
    const resultsSelect = el("select", { class: "hidden" });
    const selectedCaption = el("div", { class: "caption" });

    function refreshCaption() {
        const chosen = rankState.contexts[dbName];
        selectedCaption.textContent = chosen ? "Selected: " + chosen.context_id + " | " + chosen.context_name : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = await apiJson("/api/contexts/search?database=" + encodeURIComponent(dbName) + "&query=" + encodeURIComponent(query) + "&limit=10");
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                (data.results || []).forEach((rec, idx) => {
                    const label = rec.context_id && rec.context_name ? rec.context_id + " | " + rec.context_name : (rec.context_name || rec.context_id || "");
                    resultsSelect.appendChild(el("option", { value: String(idx) }, label));
                });
                resultsSelect._records = data.results || [];
                resultsSelect.classList.toggle("hidden", (data.results || []).length === 0);
            } catch (err) { /* keep prior dropdown contents on transient search errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        rankState.contexts[dbName] = {
            rank: rankLetter, database: dbName,
            context_id: rec.context_id || "", context_name: rec.context_name || "",
            context_type: rec.context_type || "", source: rec.source || "", reason: rec.reason || "",
        };
        renderStep3Panel(panel);
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function renderHmdbImportBlock(rankLetter, rankState, panel) {
    const wrap = el("div", {});
    wrap.appendChild(fileInputWidget({
        buttonLabel: "Choose HMDB file",
        initialFilename: rankState.hmdb ? rankState.hmdb.context_id : "",
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            formData.append("rank", rankLetter);
            try {
                const data = await apiUpload("/api/mc/hmdb-import", formData);
                rankState.hmdb = {
                    rank: rankLetter, database: "HMDB Import",
                    context_id: data.filename, context_name: data.filename, context_type: "imported_file",
                    source: "user_upload", row_count: data.row_count, key_columns: data.key_columns,
                };
                rankState.hmdbPreview = data.preview;
                renderStep3Panel(panel);
            } catch (err) {
                alert("Failed to read HMDB file: " + err.message);
            }
        },
    }));
    if (rankState.hmdb) {
        wrap.appendChild(alertBox("success", "Loaded: " + rankState.hmdb.row_count + " rows. Key columns detected: " + rankState.hmdb.key_columns.join(", ") + "."));
        wrap.appendChild(el("div", { class: "caption" }, "Preview (" + rankState.hmdb.row_count + " rows):"));
        if (rankState.hmdbPreview) wrap.appendChild(buildTable(rankState.hmdbPreview, { maxHeight: 300 }));
    }
    return wrap;
}

function renderSelectionsSummary() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "section-label" }, "Selected database/context ranks"));
    const selections = buildSelectionsList();
    STATE.mc.selections = selections;
    if (selections.length === 0) {
        wrap.appendChild(alertBox("info", "No database/context ranks selected yet."));
        return wrap;
    }
    const rows = selections.map((ctx) => ({
        "Rank": ctx.rank, "Database": ctx.database, "Context ID": ctx.context_id,
        "Context Name": ctx.context_name, "Context Type": ctx.context_type,
        "Source": ctx.database === "HMDB Import" ? "HMDB - user_upload" : ctx.database,
    }));
    wrap.appendChild(buildTable({ columns: ["Rank", "Database", "Context ID", "Context Name", "Context Type", "Source"], rows }, {}));
    return wrap;
}

function renderMcStep4() {
    const wrap = el("div", { class: "mc-results-screen" });
    wrap.appendChild(el("div", { class: "step-subheader" }, "Results"));
    const results = STATE.mc.analysisResults;
    if (!results || Object.keys(results).length === 0) {
        wrap.appendChild(alertBox("info", "Run connection analysis to generate results."));
    } else {
        Object.keys(results).forEach((rankName) => {
            const result = results[rankName];
            const rankMatch = String(rankName).match(/[A-H]/i);
            const rankLetter = rankMatch ? rankMatch[0].toLowerCase() : "a";
            const card = el("div", { class: "rank-card mc-rank-card mc-result-card mc-rank-card-" + rankLetter });
            card.appendChild(el("div", { class: "section-label mc-rank-title" }, rankName));
            if (result.selections && result.selections.length) {
                card.appendChild(buildTable({ columns: ["Database", "Context ID", "Context Name", "Context Type"], rows: result.selections }, {}));
            }
            if (result.table1 && result.table1.rows.length) {
                const det = el("details", { class: "expander" });
                det.appendChild(el("summary", {}, "Table 1 - Compound database hits"));
                const body = el("div", { class: "expander-body" });
                body.appendChild(buildTable(result.table1, {}));
                const dlBtn = el("button", { class: "btn" }, "Download table (CSV)");
                dlBtn.addEventListener("click", () => downloadFile(
                    "/api/mc/download/result?rank=" + encodeURIComponent(rankName) + "&table=table1",
                    rankName.replace(/ /g, "_") + "_compound_hits.csv",
                ));
                body.appendChild(dlBtn);
                det.appendChild(body);
                card.appendChild(det);
            }
            if (result.table2 && result.table2.rows.length) {
                const det = el("details", { class: "expander" });
                det.appendChild(el("summary", {}, "Table 2 - Rhea + UniProt enzyme presence"));
                const body = el("div", { class: "expander-body" });
                body.appendChild(buildTable(result.table2, {}));
                const dlBtn = el("button", { class: "btn" }, "Download table (CSV)");
                dlBtn.addEventListener("click", () => downloadFile(
                    "/api/mc/download/result?rank=" + encodeURIComponent(rankName) + "&table=table2",
                    rankName.replace(/ /g, "_") + "_enzyme_presence.csv",
                ));
                body.appendChild(dlBtn);
                det.appendChild(body);
                card.appendChild(det);
            }
            wrap.appendChild(card);
        });
    }

    const actions = el("div", { class: "right-actions" });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back to Rank Selection");
    backBtn.addEventListener("click", () => setWorkflowStep(3));
    const homeBtn = el("button", { class: "btn" }, "Home");
    homeBtn.addEventListener("click", () => setPage("Home"));
    btnRow.appendChild(backBtn);
    btnRow.appendChild(homeBtn);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);
    return wrap;
}

function aaHeader() {
    return el("div", { class: "aa-header" }, [
        el("img", { src: "/api/assets/logo_AA.png", alt: "Annotation Analytics" }),
        el("div", { class: "title" }, "Annotation Analytics"),
    ]);
}

function backToGenomicsActions() {
    const actions = el("div", { class: "right-actions" });
    const btn = el("button", { class: "btn" }, "Back to Genomics");
    btn.addEventListener("click", () => { STATE.aaPage = "genomics"; render(); });
    actions.appendChild(el("div", { class: "btn-row", style: { width: "25%", minWidth: "220px" } }, btn));
    return actions;
}

function renderAnnotationAnalytics() {
    if (STATE.aaPage === "genomics") return renderGenomicsAnnotationScreen();
    if (STATE.aaPage === "genomics_ora_config") return renderOraAnalysisConfigScreen();
    if (STATE.aaPage === "genomics_ora_results") return renderOraAnalysisResultsScreen();
    if (STATE.aaPage === "genomics_burden_config") return renderBurdenAnalysisConfigScreen();
    if (STATE.aaPage === "genomics_burden_results") return renderBurdenAnalysisResultsScreen();
    if (STATE.aaPage === "genomics_netprop_config") return renderNetpropAnalysisConfigScreen();
    if (STATE.aaPage === "genomics_netprop_results") return renderNetpropAnalysisResultsScreen();
    if (STATE.aaPage === "genomics_annotation_config") return renderAnnotationAnalysisConfigScreen();
    if (STATE.aaPage === "genomics_annotation_results") return renderAnnotationAnalysisResultsScreen();
    if (STATE.aaPage === "genomics_module") {
        if (STATE.genomicsModuleName === "Over-Representation Analysis") return renderGenomicsOraScreen();
        if (STATE.genomicsModuleName === "Gene-set Burden Test") return renderGenomicsBurdenScreen();
        if (STATE.genomicsModuleName === "Network Propagation / Module Detection") return renderGenomicsNetpropScreen();
        if (STATE.genomicsModuleName === "Gene Function Annotation") return renderGeneFunctionAnnotationScreen();
        return renderGenomicsUploadScreen(STATE.genomicsModuleName);
    }
    if (STATE.aaPage === "transcriptomics") return renderTranscriptomicsScreen();
    if (STATE.aaPage === "transcriptomics_module") return renderTranscriptomicsModuleScreen();
    if (STATE.aaPage === "transcriptomics_ora_config") return renderTranscriptomicsOraConfigScreen();
    if (STATE.aaPage === "transcriptomics_ora_results") return renderTranscriptomicsOraResultsScreen();
    if (STATE.aaPage === "transcriptomics_gsea_config") return renderTranscriptomicsGseaConfigScreen();
    if (STATE.aaPage === "transcriptomics_gsea_results") return renderTranscriptomicsGseaResultsScreen();
    if (STATE.aaPage === "transcriptomics_regulon_config") return renderTranscriptomicsRegulonConfigScreen();
    if (STATE.aaPage === "transcriptomics_regulon_results") return renderTranscriptomicsRegulonResultsScreen();
    if (STATE.aaPage === "transcriptomics_network_config") return renderTranscriptomicsNetworkConfigScreen();
    if (STATE.aaPage === "transcriptomics_network_results") return renderTranscriptomicsNetworkResultsScreen();
    if (STATE.aaPage === "proteomics") return renderProteomicsScreen();
    if (STATE.aaPage === "proteomics_module") return renderProteomicsModuleScreen();
    if (STATE.aaPage === "proteomics_ora_config") return renderProteomicsOraConfigScreen();
    if (STATE.aaPage === "proteomics_ora_results") return renderProteomicsOraResultsScreen();

    const wrap = el("div", {});
    wrap.appendChild(aaHeader());

    const topGrid = el("div", { class: "aa-module-grid aa-module-grid-top" });
    [["Genomics", "genomics"], ["Transcriptomics", "transcriptomics"], ["Proteomics", "proteomics"]].forEach(([name, key]) => {
        const btn = el("button", { class: "module-btn " + name.toLowerCase() }, name);
        if (key) btn.addEventListener("click", () => { STATE.aaPage = key; render(); });
        topGrid.appendChild(btn);
    });
    wrap.appendChild(topGrid);

    const bottomGrid = el("div", { class: "aa-module-grid aa-module-grid-bottom" });
    bottomGrid.appendChild(el("div", { class: "spacer" }));
    bottomGrid.appendChild(el("button", { class: "module-btn metabolomics" }, "Metabolomics"));
    bottomGrid.appendChild(el("button", { class: "module-btn multiomics" }, "Multi-Omics"));
    bottomGrid.appendChild(el("div", { class: "spacer" }));
    wrap.appendChild(bottomGrid);
    return wrap;
}

function renderGenomicsAnnotationScreen() {
    const wrap = el("div", {});
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle" }, "Genomics"));
    const grid = el("div", { class: "aa-genomics-grid" });
    [
        "Over-Representation Analysis",
        "Gene-set Burden Test",
        "Network Propagation / Module Detection",
        "Gene Function Annotation",
    ].forEach((name) => {
        const btn = el("button", { class: "module-btn genomics" }, name);
        btn.addEventListener("click", () => {
            STATE.aaPage = "genomics_module";
            STATE.genomicsModuleName = name;
            STATE.genomicsIdentifierType = "";
            STATE.ora = freshOraState();
            render();
        });
        grid.appendChild(btn);
    });
    wrap.appendChild(grid);
    return wrap;
}


// All four Transcriptomics tiles are implemented, so this set is empty - kept
// (rather than deleted) because it is the mechanism that stops a future fifth
// tile being clickable before its engine exists. A tile listed here is
// disabled up front, instead of letting a user complete a whole upload+resolve
// cycle - burning a real mygene.info job - only to find the analysis missing.
const TRANSCRIPTOMICS_TILES_NOT_YET_IMPLEMENTED = new Set([]);

function renderTranscriptomicsScreen() {
    const wrap = el("div", {});
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle" }, "Transcriptomics"));
    const grid = el("div", { class: "aa-genomics-grid" });
    [
        "Over-Representation Analysis",
        "Gene Set Enrichment Analysis",
        "Transcription Factor / Regulon",
        "Network / Module Detection",
    ].forEach((name) => {
        const notYetImplemented = TRANSCRIPTOMICS_TILES_NOT_YET_IMPLEMENTED.has(name);
        const btn = el("button", {
            class: "module-btn transcriptomics" + (notYetImplemented ? " module-btn-disabled" : ""),
            title: notYetImplemented ? name + " is not implemented yet - coming soon." : "",
        }, name + (notYetImplemented ? " (coming soon)" : ""));
        btn.disabled = notYetImplemented;
        btn.addEventListener("click", () => {
            if (notYetImplemented) return;
            STATE.aaPage = "transcriptomics_module";
            STATE.transcriptomicsModuleName = name;
            STATE.transcriptomics = freshTranscriptomicsState();
            render();
        });
        grid.appendChild(btn);
    });
    wrap.appendChild(grid);
    return wrap;
}

function backToTranscriptomicsActions() {
    const actions = el("div", { class: "right-actions" });
    const btn = el("button", { class: "btn" }, "Back to Transcriptomics");
    btn.addEventListener("click", () => { STATE.aaPage = "transcriptomics"; render(); });
    actions.appendChild(el("div", { class: "btn-row", style: { width: "25%", minWidth: "220px" } }, btn));
    return actions;
}

function renderTranscriptomicsOrganismSelector() {
    const wrap = el("div", {});
    // Network / Module Detection searches STRING's own species list instead of
    // the whole UniProt taxonomy. Same NCBI taxon-ID space, so the resolved
    // table, _require_consistent_taxon and meta["taxon"] are unaffected - but
    // STRING covers ~12,500 organisms against UniProt's ~2.4M, and picking one
    // STRING has no network for would only fail AFTER a full upload and a real
    // mygene.info resolve. The other three tiles stay unrestricted, because GO
    // and Reactome cover far more organisms than STRING does.
    const isNetworkTile = STATE.transcriptomicsModuleName === "Network / Module Detection";
    const input = el("input", { class: "text-input", type: "text",
        placeholder: isNetworkTile ? "Species / organism in STRING (name or NCBI taxon ID)"
                                   : "Species / organism (name or UniProt taxon ID)" });
    const resultsSelect = el("select", { class: "hidden" });
    const selectedCaption = el("div", { class: "caption" });

    function refreshCaption() {
        selectedCaption.textContent = STATE.transcriptomics.taxonId
            ? "Selected: " + STATE.transcriptomics.taxonName + " (taxon " + STATE.transcriptomics.taxonId + ")"
            : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = isNetworkTile
                    ? await apiJson("/api/netprop/species?query=" + encodeURIComponent(query))
                    : await apiJson("/api/contexts/search?database=" + encodeURIComponent("Rhea + UniProt") + "&query=" + encodeURIComponent(query) + "&limit=10");
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                // The two endpoints name their fields differently; normalize to
                // the {context_id, context_name} shape the change handler reads.
                const records = (data.results || []).map((rec) => (
                    isNetworkTile ? { context_id: String(rec.taxon_id || ""), context_name: String(rec.name || "") } : rec
                ));
                records.forEach((rec, idx) => {
                    const label = rec.context_id && rec.context_name ? rec.context_id + " | " + rec.context_name : (rec.context_name || rec.context_id || "");
                    resultsSelect.appendChild(el("option", { value: String(idx) }, label));
                });
                resultsSelect._records = records;
                resultsSelect.classList.toggle("hidden", records.length === 0);
            } catch (err) { /* keep prior dropdown contents on transient search errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        const newTaxonId = String(rec.context_id || "");
        // A resolved table is tied to the organism it was resolved under -
        // taxon_id drives GO/MSigDB/gene-length, taxon_name drives Reactome,
        // and switching organisms without re-resolving would leave those
        // fields (and every downstream analysis) describing two different
        // organisms at once. The server's own configure endpoints now refuse
        // that combination (_require_consistent_taxon in API.py) - this
        // clears client-side too, so the mismatch is prevented rather than
        // merely rejected after the fact.
        const hadResolvedForDifferentOrganism =
            !!STATE.transcriptomics.resolved && newTaxonId !== STATE.transcriptomics.taxonId;
        STATE.transcriptomics.taxonId = newTaxonId;
        STATE.transcriptomics.taxonName = String(rec.context_name || "");
        if (hadResolvedForDifferentOrganism) {
            _clearTranscriptomicsResolutionDependentState();
            // Shown as an inline banner (below), not a blocking alert() -
            // consistent with how every other state-changed notice in this
            // app is surfaced. Cleared once a new resolve actually succeeds
            // (renderTranscriptomicsResolveControls' onDone), not just on
            // the next render, so it survives until the user has genuinely
            // acted on it rather than flashing once and disappearing.
            STATE.transcriptomics._organismChangedNotice =
                "Organism changed to " + STATE.transcriptomics.taxonName + ". The resolved table and all " +
                "analysis results were for the previous organism, so they were cleared - re-resolve your " +
                "identifiers to continue.";
        } else {
            // Without this, picking human -> mouse (warning shown, resolved
            // cleared) -> human again (nothing left to invalidate) left the
            // FIRST switch's banner on screen, now naming the WRONG organism
            // ("changed to Mus musculus" while the current selection reads
            // Homo sapiens) - a stale notice is worse than no notice.
            STATE.transcriptomics._organismChangedNotice = null;
        }
        render();
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function renderTranscriptomicsUploadSection() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "section-label" }, "Upload differential expression table"));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.transcriptomics.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            formData.append("identifier_type", STATE.transcriptomics.identifierType);
            try {
                const data = await apiUpload("/api/transcriptomics/upload", formData);
                STATE.transcriptomics.deTable = data.table;
                STATE.transcriptomics.uploadedFilename = file.name;
                STATE.transcriptomics.columnMapping = data.column_mapping;
                // The SERVER's upload endpoint replaces the entire session
                // sub-dict (_fresh_transcriptomics_state() in API.py), so
                // ora/gsea/regulon config+results are already gone server-
                // side the moment this response arrives - clear them here
                // too, or the client would keep showing results the server
                // can no longer serve (e.g. a CSV download would 400).
                STATE.transcriptomics.customDb = null;
                STATE.transcriptomics.background = null;
                _clearTranscriptomicsResolutionDependentState();
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));

    const idLabel = TRANSCRIPTOMICS_ID_LABELS[STATE.transcriptomics.identifierType] || "identifier";
    wrap.appendChild(el("div", { class: "caption" },
        "Required columns: " + idLabel + " + log2FoldChange + pvalue. Optional column: padj/FDR " +
        "(computed automatically via Benjamini-Hochberg if omitted)."));

    const exLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Example file");
    exLink.addEventListener("click", () => {
        const idType = STATE.transcriptomics.identifierType;
        downloadFile(
            "/api/transcriptomics/example?identifier_type=" + encodeURIComponent(idType),
            TRANSCRIPTOMICS_EXAMPLE_FILENAMES[idType] || "example_transcriptomics.csv",
        );
    });
    wrap.appendChild(exLink);
    // The transcript example is intentionally tiny - it demonstrates
    // transcript COLLAPSING, which is a retriever behaviour, so enlarging it
    // would defeat its purpose. Say so, rather than leaving a 7-gene example
    // looking like an oversight next to the ~3000-gene ones.
    if (STATE.transcriptomics.identifierType === "ensembl_transcript") {
        wrap.appendChild(el("div", { class: "caption" },
            "The Ensembl Transcript example is deliberately small - it exists to show how multiple " +
            "transcripts of the same gene are collapsed into one row. It is too small for the " +
            "downstream analyses to return meaningful results; use the Gene Symbol, Ensembl Gene ID " +
            "or Entrez Gene ID example for those."));
    }

    if (STATE.transcriptomics.deTable) {
        if (STATE.transcriptomics.columnMapping) {
            const cm = STATE.transcriptomics.columnMapping;
            const parts = Object.keys(cm).map((k) => k + " -> " + cm[k]);
            wrap.appendChild(el("div", { class: "caption" }, "Detected columns: " + parts.join(", ")));
        }
        const maxHeight = Math.min(400, Math.max(200, (STATE.transcriptomics.deTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.transcriptomics.deTable, { maxHeight }));
    }

    return wrap;
}

function renderTranscriptomicsResolveControls() {
    const wrap = el("div", {});
    const disabled = !STATE.transcriptomics.deTable || !STATE.transcriptomics.taxonId;
    const btn = el("button", { class: "btn" }, "Resolve identifiers");
    btn.disabled = disabled;
    btn.addEventListener("click", async () => {
        btn.disabled = true;
        const spinner = el("div", { class: "caption" }, "Resolving identifiers via mygene.info...");
        wrap.appendChild(spinner);
        try {
            const { job_id } = await apiPostJson("/api/transcriptomics/resolve", {
                taxon_id: STATE.transcriptomics.taxonId,
                taxon_name: STATE.transcriptomics.taxonName,
            });
            pollJob(job_id, {
                onDone: (result) => {
                    STATE.transcriptomics.resolved = result;
                    STATE.transcriptomics._organismChangedNotice = null;
                    render();
                },
                onError: (err) => { alert(err); render(); },
            });
        } catch (err) {
            alert(err.message);
            render();
        }
    });
    wrap.appendChild(btn);
    return wrap;
}

function renderTranscriptomicsResults() {
    const wrap = el("div", {});
    const resolved = STATE.transcriptomics.resolved;
    const meta = resolved.metrics;
    wrap.appendChild(el("div", { class: "section-label" }, "Resolved differential expression table"));
    wrap.appendChild(buildTable(resolved.resolved_table, { maxHeight: 420 }));
    wrap.appendChild(metricCards([
        ["Mapped", meta.n_mapped],
        ["Unmapped", meta.n_unmapped],
        ["Duplicates collapsed", meta.n_duplicates_collapsed],
        ["Transcripts collapsed", meta.n_transcripts_collapsed],
        ["Dropped (non-numeric)", meta.n_dropped_nonnumeric],
    ], 5));
    if (meta.padj_computed) {
        wrap.appendChild(alertBox("info",
            "No adjusted p-value (padj/FDR) column was found - it was computed automatically via " +
            "Benjamini-Hochberg from the p-value column."));
    }
    if (meta.n_recovered_fallback > 0) {
        wrap.appendChild(alertBox("info",
            meta.n_recovered_fallback + " gene(s) recovered by a broadened identifier search " +
            "(e.g. deprecated/retired gene symbols) after the primary lookup found no match."));
    }
    if (meta.n_ambiguous > 0) {
        wrap.appendChild(alertBox("info",
            meta.n_ambiguous + " identifier(s) matched more than one gene; the first match " +
            "returned was kept for each."));
    }
    if (meta.n_fetch_failures > 0) {
        wrap.appendChild(alertBox("warning",
            meta.n_fetch_failures + " mygene.info request(s) failed after retries" +
            (meta.n_ids_affected_by_failure > 0
                ? ", leaving " + meta.n_ids_affected_by_failure + " identifier(s) unresolved for reasons " +
                  "unrelated to whether those genes exist (a network issue, not \"not found\")"
                : "") +
            ". Consider re-running if this table looks incomplete."));
    }

    const dlResolvedBtn = el("button", { class: "btn" }, "Download resolved table (CSV)");
    dlResolvedBtn.addEventListener("click", () => {
        downloadFile("/api/transcriptomics/download/resolved.csv", "transcriptomics_resolved.csv");
    });
    wrap.appendChild(dlResolvedBtn);

    if (STATE.transcriptomicsModuleName === "Over-Representation Analysis") {
        const continueBtn = el("button", { class: "btn btn-primary", style: { marginTop: "0.8rem" } }, "Continue to Over-Representation Analysis");
        continueBtn.addEventListener("click", () => {
            STATE.transcriptomics.oraConfig = null;
            STATE.transcriptomics.oraResults = null;
            STATE.transcriptomics.oraRunLocked = false;
            STATE.aaPage = "transcriptomics_ora_config";
            render();
        });
        wrap.appendChild(continueBtn);
    } else if (STATE.transcriptomicsModuleName === "Gene Set Enrichment Analysis") {
        const continueBtn = el("button", { class: "btn btn-primary", style: { marginTop: "0.8rem" } }, "Continue to Gene Set Enrichment Analysis");
        continueBtn.addEventListener("click", () => {
            STATE.transcriptomics.gseaConfig = null;
            STATE.transcriptomics.gseaResults = null;
            STATE.transcriptomics.gseaRunLocked = false;
            STATE.aaPage = "transcriptomics_gsea_config";
            render();
        });
        wrap.appendChild(continueBtn);
    } else if (STATE.transcriptomicsModuleName === "Transcription Factor / Regulon") {
        const continueBtn = el("button", { class: "btn btn-primary", style: { marginTop: "0.8rem" } }, "Continue to Transcription Factor / Regulon");
        continueBtn.addEventListener("click", () => {
            STATE.transcriptomics.regulonConfig = null;
            STATE.transcriptomics.regulonResults = null;
            STATE.transcriptomics.regulonRunLocked = false;
            STATE.aaPage = "transcriptomics_regulon_config";
            render();
        });
        wrap.appendChild(continueBtn);
    } else if (STATE.transcriptomicsModuleName === "Network / Module Detection") {
        const continueBtn = el("button", { class: "btn btn-primary", style: { marginTop: "0.8rem" } }, "Continue to Network / Module Detection");
        continueBtn.addEventListener("click", () => {
            STATE.transcriptomics.networkConfig = null;
            STATE.transcriptomics.networkResults = null;
            STATE.transcriptomics.networkRunLocked = false;
            STATE.aaPage = "transcriptomics_network_config";
            render();
        });
        wrap.appendChild(continueBtn);
    } else {
        // Defensive: every tile in the grid above now has a branch, so this is
        // only reachable if a new tile is added without one.
        wrap.appendChild(el("div", { class: "caption" },
            "This analysis module is not available yet."));
    }
    return wrap;
}

function renderTranscriptomicsModuleScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — " + STATE.transcriptomicsModuleName));

    const typeRow = el("div", { class: "row" });
    function typeBtn(label, type) {
        const btn = el("button", { class: "btn" }, label);
        if (STATE.transcriptomics.identifierType === type) btn.style.fontWeight = "800";
        btn.addEventListener("click", () => {
            if (STATE.transcriptomics.identifierType !== type) {
                const moduleName = STATE.transcriptomicsModuleName;
                STATE.transcriptomics = freshTranscriptomicsState();
                STATE.transcriptomicsModuleName = moduleName;
                STATE.transcriptomics.identifierType = type;
                render();
            }
        });
        return btn;
    }
    [["Gene Symbol", "gene_symbol"], ["Ensembl Gene ID", "ensembl_gene"],
     ["Ensembl Transcript ID", "ensembl_transcript"], ["Entrez Gene ID", "entrez"]].forEach(([label, type]) => {
        typeRow.appendChild(el("div", { class: "col" }, typeBtn(label, type)));
    });
    wrap.appendChild(typeRow);

    if (!STATE.transcriptomics.identifierType) {
        wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
        wrap.appendChild(backToTranscriptomicsActions());
        return wrap;
    }

    wrap.appendChild(el("div", { class: "section-label" }, "Select organism"));
    wrap.appendChild(renderTranscriptomicsOrganismSelector());
    if (STATE.transcriptomics._organismChangedNotice) {
        wrap.appendChild(alertBox("warning", STATE.transcriptomics._organismChangedNotice));
    }

    wrap.appendChild(renderTranscriptomicsUploadSection());

    wrap.appendChild(renderTranscriptomicsResolveControls());

    if (STATE.transcriptomics.resolved) {
        wrap.appendChild(renderTranscriptomicsResults());
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// Transcriptomics — Over-Representation Analysis (ORA)
// ---------------------------------------------------------------------------

// Collections offered from MSigDB (CC BY 4.0, (c) Broad Institute/MIT/UC
// Regents - see info/README-ORA-transcriptomics.txt "MSigDB collections and
// licensing" for the full list, what was excluded, and why). Keys and order
// here MUST match Transcriptomics_module.MSIGDB_COLLECTIONS exactly - a
// Playwright check guards this. mouseAvailable=false means the collection
// genuinely has no MSigDB file for mouse (verified against MSigDB's own
// collections page), not merely unfetched.
const TRANSCRIPTOMICS_MSIGDB_COLLECTIONS = [
    { key: "hallmark", label: "Hallmark", mouseAvailable: true },
    { key: "perturbations", label: "Perturbations (CGP)", mouseAvailable: true },
    { key: "wikipathways", label: "WikiPathways", mouseAvailable: true },
    { key: "pid", label: "PID", mouseAvailable: false },
    { key: "phenotype", label: "Phenotype (HPO / MPT)", mouseAvailable: true },
    { key: "oncogenic", label: "Oncogenic", mouseAvailable: false },
    { key: "immunologic", label: "Immunologic", mouseAvailable: true },
    { key: "cell_type", label: "Cell type", mouseAvailable: true },
];
// MSigDB covers human and mouse only - a different, wider organism whitelist
// than the per-collection mouseAvailable flags above, which only matter once
// the organism is confirmed to be one of these two.
function _msigdbAvailableForTaxon(taxonId) {
    return taxonId === "9606" || taxonId === "10090";
}
function _msigdbCollectionAvailableForTaxon(collection, taxonId) {
    if (taxonId === "9606") return true;
    if (taxonId === "10090") return collection.mouseAvailable;
    return false;
}

// Shown only when the "MSigDB" pill is checked - lets the user pick which of
// the 8 offered collections to run. EACH selected collection becomes its own
// independently BH-corrected result family (see the family-count caption on
// the results screen) - shared between the ORA and GSEA config screens.
function renderMsigdbCollectionSelector(locked, draft, taxonId, onChange) {
    const wrap = el("div", { class: "ora-msigdb-collections" });
    wrap.appendChild(el("div", { class: "chart-option-label" }, "MSigDB collections"));
    const row = el("div", { class: "db-checkbox-row" });
    let anyUnavailable = false;
    TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.forEach((c) => {
        const available = _msigdbCollectionAvailableForTaxon(c, taxonId);
        if (!available) {
            anyUnavailable = true;
            draft.msigdbCollections = draft.msigdbCollections.filter((k) => k !== c.key);
        }
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked || !available;
        cb.checked = available && draft.msigdbCollections.indexOf(c.key) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.msigdbCollections.indexOf(c.key) === -1) draft.msigdbCollections.push(c.key);
            } else {
                draft.msigdbCollections = draft.msigdbCollections.filter((k) => k !== c.key);
            }
            onChange();
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + c.label + (available ? "" : " (no mouse data)")));
        row.appendChild(label);
    });
    wrap.appendChild(row);
    if (anyUnavailable) {
        wrap.appendChild(el("div", { class: "caption" },
            "Some collections have no MSigDB file for mouse at all (not merely unfetched) - " +
            "see info/README-ORA-transcriptomics.txt for the full per-collection availability table."));
    }
    if (!draft.msigdbCollections.length) {
        wrap.appendChild(el("div", { class: "caption" }, "Select at least one collection to run MSigDB."));
    }
    return wrap;
}

// ---------------------------------------------------------------------------
// MSigDB cross-collection summary (shared core, ORA + GSEA-specific wrappers)
// - the ONE genuinely new figure this feature needs. With up to 8 collection
// families per direction, nothing else answers "which collection is
// carrying the signal?" without clicking through every rail entry by hand.
// Every other MSigDB result reuses the dot/bar/volcano/enrichment-map/
// gene-concept-network builders every other database already has.
// ---------------------------------------------------------------------------
function _buildMsigdbCollectionSummaryCore(div, families, alpha, slug, opts) {
    if (typeof Plotly === "undefined") return;
    if (!families.length) { div.textContent = "No MSigDB collections in this run."; return; }
    opts = opts || {};
    // GSEA's NES is genuinely SIGNED (positive = enriched toward the up end,
    // negative toward the down end), so colouring by its sign is meaningful.
    // ORA's is NOT: its test is one-sided over-representation only, its effect
    // measure is fold enrichment (a non-negative ratio whose neutral point is
    // 1, not 0), and the "(all)" family pools up- and down-regulated DEGs
    // anyway. An earlier version of this figure applied `bestEffect >= 0` to
    // BOTH - which for ORA is true for every possible value, so the blue
    // branch was unreachable and the colour silently encoded nothing while
    // looking like it encoded direction. ORA now uses one enrichment colour.
    const directional = !!opts.directional;
    const effectLabel = opts.effectLabel || "effect";
    // Canonical order (TRANSCRIPTOMICS_MSIGDB_COLLECTIONS), reversed so the
    // first collection lands at the TOP of the horizontal bar (Plotly draws
    // categorical y-axes bottom-up by default).
    families = families.slice().reverse();
    const labels = families.map((f) => f.label);
    const counts = families.map((f) => f.nSignificant);
    const colors = families.map((f) => {
        if (f.nSignificant === 0) return "#b8c2cc";
        if (!directional) return "#c0392b";
        // An unknown best effect is NOT silently painted "down" - grey it,
        // the same as "nothing significant", since its direction is unknown.
        if (f.bestEffect === null || typeof f.bestEffect !== "number" || !isFinite(f.bestEffect)) return "#b8c2cc";
        return f.bestEffect >= 0 ? "#c0392b" : "#2874a6";
    });
    const _fin = (v) => v !== null && typeof v === "number" && isFinite(v);
    const fmt = (v, digits) => _fin(v) ? v.toFixed(digits) : "n/a";
    // See _fmtSci's own comment for why "0" is shown as "< 1e-300" rather
    // than "0.00e+0". _fin() already excludes null/non-finite, so the only
    // difference from _fmtSci here is the "n/a" fallback text.
    const fmtExp = (v) => _fin(v) ? _fmtSci(v, 2) : "n/a";
    const trace = {
        x: counts, y: labels, type: "bar", orientation: "h",
        marker: { color: colors },
        customdata: families.map((f) => [
            fmtExp(f.bestFdr),
            fmt(f.bestEffect, 2),
            f.bestTermName || "n/a",
        ]),
        hovertemplate: "<b>%{y}</b><br>Significant sets: %{x}<br>Best FDR: %{customdata[0]}" +
            "<br>Best " + effectLabel + ": %{customdata[1]}" +
            "<br>Top term: %{customdata[2]}<extra></extra>",
    };
    const layout = {
        xaxis: { title: { text: "Significant gene sets (FDR < " + alpha + ")" }, zeroline: true, rangemode: "tozero" },
        yaxis: { automargin: true },
        height: Math.max(220, families.length * 42 + 60),
        margin: { l: 170, r: 30, t: 20, b: 46 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "msigdbsummary", "MSigDB collection summary");
}

// families: [{label, nSignificant, bestFdr, bestEffect, bestTermName}, ...]
// extracted from result tables, in TRANSCRIPTOMICS_MSIGDB_COLLECTIONS' own
// canonical order - the extraction (which columns/keys to read) is the only
// part that differs between ORA and GSEA; the rendering above is identical.
// Strict numeric coercion for result-table cells. Plain Number() is WRONG
// here: Number(null), Number("") and Number(" ") are all 0, not NaN, and a
// pandas NaN arrives over JSON as null - so an unguarded Number() turns a
// MISSING FDR into 0, i.e. silently the most significant value possible,
// which would then be reported as the family's top hit and counted as
// significant. Only genuine finite numbers get through here.
function _msigdbNum(v) {
    if (v === null || v === undefined || v === "") return NaN;
    const n = Number(v);
    return isFinite(n) ? n : NaN;
}

// Shared "best row" scan. The NaN guard is load-bearing, not defensive
// noise: `NaN < bestFdr` is false for every comparison, so seeding bestFdr
// from an unguarded first row would pin it to NaN forever and report that
// row as the family's top hit no matter what followed it.
function _msigdbBestRow(rows, fdrKey, effectKey) {
    let bestFdr = null, bestEffect = null, bestTermName = null;
    rows.forEach((r) => {
        const fdr = _msigdbNum(r[fdrKey]);
        if (isNaN(fdr)) return;
        if (bestFdr === null || fdr < bestFdr) {
            bestFdr = fdr;
            const eff = _msigdbNum(r[effectKey]);
            bestEffect = isNaN(eff) ? null : eff;
            bestTermName = r["Term Name"];
        }
    });
    return { bestFdr, bestEffect, bestTermName };
}

function _msigdbFamiliesFromOraTables(tables, alpha) {
    const families = [];
    TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.forEach((c) => {
        const key = "MSigDB: " + c.label + " (all)";
        const block = tables[key];
        if (!block) return;
        const rows = _tableToObjects(block.table);
        const nSignificant = rows.filter((r) => {
            const fdr = _msigdbNum(r.FDR);
            return !isNaN(fdr) && fdr < alpha;
        }).length;
        const best = _msigdbBestRow(rows, "FDR", "Fold Enrichment");
        families.push({ label: c.label, nSignificant, ...best });
    });
    return families;
}
function _msigdbFamiliesFromGseaTables(tables, alpha) {
    const families = [];
    TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.forEach((c) => {
        const key = "MSigDB: " + c.label;
        const block = tables[key];
        if (!block) return;
        const rows = _tableToObjects(block.table);
        const nSignificant = rows.filter((r) => {
            const fdr = _msigdbNum(r["FDR q-value"]);
            return !isNaN(fdr) && fdr < alpha;
        }).length;
        const best = _msigdbBestRow(rows, "FDR q-value", "NES");
        families.push({ label: c.label, nSignificant, ...best });
    });
    return families;
}
function buildMsigdbOraCollectionSummary(div, tables, alpha, slug) {
    // directional:false - ORA's fold enrichment is unsigned (see the note in
    // _buildMsigdbCollectionSummaryCore).
    _buildMsigdbCollectionSummaryCore(div, _msigdbFamiliesFromOraTables(tables, alpha), alpha, slug,
        { directional: false, effectLabel: "fold enrichment" });
}
function buildMsigdbGseaCollectionSummary(div, tables, alpha, slug) {
    _buildMsigdbCollectionSummaryCore(div, _msigdbFamiliesFromGseaTables(tables, alpha), alpha, slug,
        { directional: true, effectLabel: "NES" });
}

const TRANSCRIPTOMICS_ORA_DATABASES = ["GO", "Reactome", "MSigDB", "Custom database"];

function renderTranscriptomicsOraConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Over-Representation Analysis: database"));

    const locked = STATE.transcriptomics.oraRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    if (!STATE.transcriptomics._oraDraft) {
        STATE.transcriptomics._oraDraft = {
            databases: [],
            msigdbCollections: [],
            padjCutoff: 0.05,
            lfcCutoff: 1.0,
            // clusterProfiler's defaults (min 10 / max 500) - matches
            // Transcriptomics_module.py's MIN/MAX_TERM_SIZE_DEFAULT.
            minTermSize: 10,
            maxTermSize: 500,
            useLengthCorrection: true,
            includeIea: true,
            useBackground: false,
        };
    }
    const draft = STATE.transcriptomics._oraDraft;

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const runBtnWrap = el("div", { class: "col" });

    function renderOraRunBtn() {
        runBtnWrap.innerHTML = "";
        const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
        const needsCustomDb = draft.databases.indexOf("Custom database") !== -1 && !STATE.transcriptomics.customDb;
        const needsMsigdbCollection = draft.databases.indexOf("MSigDB") !== -1 && draft.msigdbCollections.length === 0;
        runBtn.disabled = locked || draft.databases.length === 0 || needsCustomDb || needsMsigdbCollection;
        runBtn.addEventListener("click", () => startTranscriptomicsOraRun(draft));
        runBtnWrap.appendChild(runBtn);
    }

    // MSigDB covers human and mouse only - for any other organism the
    // checkbox is disabled (and unchecked) so a run isn't misled by
    // silently-empty results.
    const taxonId = STATE.transcriptomics.taxonId;
    const msigdbAvailable = _msigdbAvailableForTaxon(taxonId);
    TRANSCRIPTOMICS_ORA_DATABASES.forEach((dbName) => {
        const unavailable = dbName === "MSigDB" && !msigdbAvailable;
        if (unavailable) {
            draft.databases = draft.databases.filter((d) => d !== dbName);
        }
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked || unavailable;
        cb.checked = !unavailable && draft.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.databases.indexOf(dbName) === -1) draft.databases.push(dbName);
            } else {
                draft.databases = draft.databases.filter((d) => d !== dbName);
            }
            // "Custom database" toggles a whole upload widget below, and
            // "MSigDB" toggles the collection sub-selector, not just the run
            // button, so both need a full re-render; the others don't.
            if (dbName === "Custom database" || dbName === "MSigDB") render();
            else renderOraRunBtn();
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName + (unavailable ? " (human/mouse only)" : "")));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    if (!msigdbAvailable) {
        wrap.appendChild(el("div", { class: "caption" },
            "MSigDB is disabled: it covers human and mouse only. GO and Reactome support this organism."));
    }
    if (draft.databases.indexOf("MSigDB") !== -1 && msigdbAvailable) {
        wrap.appendChild(renderMsigdbCollectionSelector(locked, draft, taxonId, renderOraRunBtn));
    }
    if (draft.databases.indexOf("Custom database") !== -1) {
        wrap.appendChild(renderTranscriptomicsCustomDbWidget(locked));
    }
    wrap.appendChild(el("div", { class: "caption" },
        "GO: Biological Process / Molecular Function / Cellular Component (EBI QuickGO, tested separately). " +
        "Reactome: pathway gene sets (native Reactome file). MSigDB: pick from 8 curated collections " +
        "(CC BY 4.0), each its own result family. Custom database: your own uploaded gene sets. Each is " +
        "tested against all DEGs, and separately against up- and down-regulated DEGs."));

    wrap.appendChild(renderTranscriptomicsBackgroundWidget(locked, draft));

    function numberField(labelText, value, step, onChange) {
        const col = el("div", { class: "col" });
        col.appendChild(el("div", { class: "chart-option-label" }, labelText));
        const input = el("input", { type: "number", class: "text-input", step: step });
        input.value = value;
        input.disabled = locked;
        input.addEventListener("change", () => {
            const v = parseFloat(input.value);
            if (isFinite(v)) onChange(v);
        });
        col.appendChild(input);
        return col;
    }

    const thresholdRow = el("div", { class: "row", style: { marginTop: "0.8rem" } });
    thresholdRow.appendChild(numberField("Adjusted p-value (padj) cutoff", draft.padjCutoff, "0.001", (v) => { draft.padjCutoff = v; }));
    thresholdRow.appendChild(numberField("|log2 Fold-Change| cutoff", draft.lfcCutoff, "0.1", (v) => { draft.lfcCutoff = v; }));
    wrap.appendChild(thresholdRow);

    const termSizeRow = el("div", { class: "row" });
    const minCol = el("div", { class: "col" });
    minCol.appendChild(el("div", { class: "chart-option-label" }, "Minimum term size (genes)"));
    const minInput = el("input", { type: "number", class: "text-input", min: "1", step: "1" });
    minInput.value = draft.minTermSize;
    minInput.disabled = locked;
    minInput.addEventListener("change", () => {
        const v = parseInt(minInput.value, 10);
        draft.minTermSize = (isFinite(v) && v >= 1) ? v : 10;
        minInput.value = draft.minTermSize;
    });
    minCol.appendChild(minInput);
    termSizeRow.appendChild(minCol);

    const maxCol = el("div", { class: "col" });
    maxCol.appendChild(el("div", { class: "chart-option-label" }, "Maximum term size (genes, blank = no limit)"));
    const maxInput = el("input", { type: "number", class: "text-input", min: "1", step: "1", placeholder: "no limit" });
    maxInput.value = draft.maxTermSize;
    maxInput.disabled = locked;
    maxInput.addEventListener("change", () => {
        const raw = maxInput.value.trim();
        draft.maxTermSize = raw === "" ? "" : parseInt(raw, 10);
    });
    maxCol.appendChild(maxInput);
    termSizeRow.appendChild(maxCol);
    wrap.appendChild(termSizeRow);

    const lcRow = el("div", { class: "row", style: { alignItems: "center", marginTop: "0.5rem" } });
    const lcLabel = el("label", { style: { display: "flex", alignItems: "center", gap: "0.4rem" } });
    const lcCheckbox = el("input", { type: "checkbox" });
    lcCheckbox.checked = draft.useLengthCorrection;
    lcCheckbox.disabled = locked;
    lcCheckbox.addEventListener("change", () => { draft.useLengthCorrection = lcCheckbox.checked; });
    lcLabel.appendChild(lcCheckbox);
    lcLabel.appendChild(document.createTextNode(" Correct for gene-length bias (Wallenius/GOseq-style)"));
    lcRow.appendChild(lcLabel);
    wrap.appendChild(lcRow);
    wrap.appendChild(el("div", { class: "caption" },
        "RNA-seq detects longer genes more easily, independent of biology. When enabled, terms are tested " +
        "with a length-bias-corrected statistic; when gene-length data is unavailable for the selected " +
        "organism, the analysis automatically falls back to the standard (uncorrected) test, and this is " +
        "reported after the run."));

    const ieaRow = el("div", { class: "row", style: { alignItems: "center", marginTop: "0.5rem" } });
    const ieaLabel = el("label", { style: { display: "flex", alignItems: "center", gap: "0.4rem" } });
    const ieaCheckbox = el("input", { type: "checkbox", class: "transcriptomics-ora-iea-checkbox" });
    ieaCheckbox.checked = draft.includeIea;
    ieaCheckbox.disabled = locked;
    ieaCheckbox.addEventListener("change", () => { draft.includeIea = ieaCheckbox.checked; });
    ieaLabel.appendChild(ieaCheckbox);
    ieaLabel.appendChild(document.createTextNode(" Include electronic annotations (IEA)"));
    ieaRow.appendChild(ieaLabel);
    wrap.appendChild(ieaRow);
    wrap.appendChild(el("div", { class: "caption" },
        "IEA (\"Inferred from Electronic Annotation\") terms are assigned automatically, without curator " +
        "review - most GO annotations are IEA, and Reactome carries the same evidence distinction for its " +
        "orthology-projected pathways. Uncheck to restrict both GO and Reactome to curated evidence only. " +
        "Has no effect on MSigDB, which carries no evidence-code data at all. Note: Reactome's " +
        "curated (non-IEA) coverage is essentially human-only - unchecking this for a non-human organism " +
        "will likely leave Reactome with zero pathways."));

    const actions = el("div", { class: "right-actions", style: { marginTop: "1rem" } });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "transcriptomics_module"; render(); });
    btnRow.appendChild(backBtn);
    renderOraRunBtn();
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function renderTranscriptomicsCustomDbWidget(locked) {
    const wrap = el("div", { class: "ora-custom-db-widget" });
    wrap.appendChild(el("div", { class: "section-label" }, "Custom database"));
    wrap.appendChild(el("div", { class: "caption" },
        "Upload your own gene sets (GMT or a tidy term/gene table), always in gene-symbol space."));

    if (!locked) {
        wrap.appendChild(fileInputWidget({
            accept: ".gmt,.txt,.tsv,.csv",
            buttonLabel: "Choose gene-set file",
            initialFilename: STATE.transcriptomics.customDb ? STATE.transcriptomics.customDb.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                try {
                    const data = await apiUpload("/api/transcriptomics/ora/custom-db-upload", formData);
                    STATE.transcriptomics.customDb = data;
                    render();
                } catch (err) {
                    STATE.transcriptomics.customDb = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    if (STATE.transcriptomics.customDb) {
        const s = STATE.transcriptomics.customDb;
        wrap.appendChild(el("div", { class: "caption ora-custom-summary" },
            "Loaded \"" + s.filename + "\": " + s.n_terms + " gene set(s), " + s.universe_size + " unique genes."));
    }

    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "(example input file)");
    exampleLink.addEventListener("click", () => {
        downloadFile("/api/transcriptomics/ora/custom-db-example", "example_custom_geneset_db_symbols.gmt");
    });
    wrap.appendChild(el("div", { class: "caption ora-custom-help" }, [
        "Format: GMT (Term | Description | Gene1 ...) or a tidy table with columns term_id, term_name, gene. " +
        "The background universe is the union of every gene in the file, intersected with the genes tested " +
        "in this experiment (same convention as GO/Reactome/Hallmark above). ",
        exampleLink,
    ]));

    return wrap;
}

// Optional background (reference) gene list, applied globally to ALL
// selected databases (unlike Genomics, where background is per-database) -
// this module already has one shared universe (tested genes), so there is
// one shared background restriction on top of it. Checkbox-gated, same
// collapsed-by-default pattern as Genomics' renderOraBackgroundWidget.
function renderTranscriptomicsBackgroundWidget(locked, draft) {
    const wrap = el("div", { class: "ora-background-widget" });

    const toggleRow = el("div", { class: "ora-iea-row ora-background-toggle-row" });
    const toggleCb = el("input", { type: "checkbox", class: "ora-iea-checkbox ora-background-checkbox" });
    toggleCb.checked = draft.useBackground;
    toggleCb.disabled = locked;
    toggleCb.addEventListener("change", () => {
        draft.useBackground = toggleCb.checked;
        render();
    });
    toggleRow.appendChild(toggleCb);
    toggleRow.appendChild(document.createTextNode(" Background gene list (optional)"));
    wrap.appendChild(toggleRow);

    if (!draft.useBackground) return wrap;

    const body = el("div", { class: "ora-background-body" });

    if (!locked) {
        body.appendChild(fileInputWidget({
            accept: ".txt,.csv,.tsv",
            buttonLabel: "Choose background file",
            initialFilename: STATE.transcriptomics.background ? STATE.transcriptomics.background.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                try {
                    const data = await apiUpload("/api/transcriptomics/ora/background-upload", formData);
                    STATE.transcriptomics.background = { filename: data.filename, n_genes: data.n_genes };
                    render();
                } catch (err) {
                    STATE.transcriptomics.background = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    if (STATE.transcriptomics.background) {
        const s = STATE.transcriptomics.background;
        body.appendChild(el("div", { class: "caption ora-background-summary" },
            "Background: \"" + s.filename + "\" - " + s.n_genes + " gene symbols. " +
            "Restricts the universe (for every selected database) to genes present in both this list and the " +
            "genes tested in this experiment."));
        if (!locked) {
            const removeBtn = el("button", { class: "btn btn-small ora-background-remove", type: "button" }, "Remove background");
            removeBtn.addEventListener("click", async () => {
                try { await apiPostJson("/api/transcriptomics/ora/background-remove", {}); } catch (err) { /* ignore */ }
                STATE.transcriptomics.background = null;
                render();
            });
            body.appendChild(removeBtn);
        }
    } else {
        const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "example file");
        exampleLink.addEventListener("click", () => {
            downloadFile("/api/transcriptomics/ora/background-example", "example_transcriptomics_background.csv");
        });
        body.appendChild(el("div", { class: "caption ora-background-help" }, [
            "Provide gene symbols (one per line, or a column headed \"Gene Symbol\") to use a custom reference set " +
            "(e.g. every gene expressed in your assay). The universe becomes the genes present in both this list and " +
            "your results table. Leave this off to use the default: every tested gene. (",
            exampleLink,
            ")",
        ]));
    }

    wrap.appendChild(body);
    return wrap;
}

function startTranscriptomicsOraRun(draft) {
    if (draft.databases.length === 0) return;
    const payload = {
        // Only sent when MSigDB is actually selected. The draft deliberately
        // REMEMBERS collection picks across an uncheck/re-check of the MSigDB
        // pill (convenient), so without this guard a run with MSigDB
        // unselected would still post collection keys the backend then stores
        // in its echoed config while never running them - config that claims
        // more than the run did.
        msigdb_collections: draft.databases.indexOf("MSigDB") !== -1 ? draft.msigdbCollections : [],
        databases: draft.databases,
        padj_cutoff: draft.padjCutoff,
        lfc_cutoff: draft.lfcCutoff,
        min_term_size: draft.minTermSize,
        max_term_size: draft.maxTermSize === "" ? null : draft.maxTermSize,
        use_length_correction: draft.useLengthCorrection,
        include_iea: draft.includeIea,
        use_background: !!(draft.useBackground && STATE.transcriptomics.background),
    };
    apiPostJson("/api/transcriptomics/ora/configure", payload).then(() => {
        STATE.transcriptomics.oraConfig = payload;
        STATE.transcriptomics.oraRunLocked = true;
        STATE.transcriptomics.oraResults = null;
        STATE.aaPage = "transcriptomics_ora_results";
        render();
        return apiPostJson("/api/transcriptomics/ora/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("transcriptomics-ora-progress-inner");
                const label = document.getElementById("transcriptomics-ora-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => {
                STATE.transcriptomics.oraResults = result;
                render();
            },
            onError: (err) => {
                alert("ORA run failed: " + err);
                STATE.transcriptomics.oraRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.transcriptomics.oraRunLocked = false;
        render();
    });
}

const TRANSCRIPTOMICS_ORA_FAMILY_ORDER = [
    "GO: Biological Process", "GO: Molecular Function", "GO: Cellular Component",
    "Reactome",
    ...TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.map((c) => "MSigDB: " + c.label),
    "Custom database",
];
const TRANSCRIPTOMICS_ORA_DIRECTION_ORDER = ["all", "up", "down"];

function _splitTranscriptomicsOraKey(key) {
    const m = key.match(/^(.*) \((all|up|down)\)$/);
    return m ? [m[1], m[2]] : [key, "all"];
}

function _sortedTranscriptomicsOraKeys(tables) {
    return Object.keys(tables).sort((a, b) => {
        const [famA, dirA] = _splitTranscriptomicsOraKey(a);
        const [famB, dirB] = _splitTranscriptomicsOraKey(b);
        const iA = TRANSCRIPTOMICS_ORA_FAMILY_ORDER.indexOf(famA);
        const iB = TRANSCRIPTOMICS_ORA_FAMILY_ORDER.indexOf(famB);
        if (iA !== iB) return iA - iB;
        return TRANSCRIPTOMICS_ORA_DIRECTION_ORDER.indexOf(dirA) - TRANSCRIPTOMICS_ORA_DIRECTION_ORDER.indexOf(dirB);
    });
}

// Result keys look like "GO: Biological Process (all)" - unsafe as a chart
// element id or an export filename (":" is invalid on Windows). Slugify once
// per key and pass the slug to every chart builder instead of the raw key.
function _transcriptomicsOraSlug(key) {
    return key.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "");
}

function renderTranscriptomicsOraResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Over-Representation Analysis: results"));

    if (!STATE.transcriptomics.oraResults) {
        wrap.appendChild(alertBox("info", "Running analysis..."));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "transcriptomics-ora-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "transcriptomics-ora-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const { tables, summary } = STATE.transcriptomics.oraResults;

    wrap.appendChild(metricCards([
        ["Universe (tested genes)", summary.n_universe],
        ["DEGs (all)", summary.n_deg_all],
        ["DEGs (up)", summary.n_deg_up],
        ["DEGs (down)", summary.n_deg_down],
    ], 4));

    if (summary.background_applied) {
        const degsDropped = summary.n_deg_all_before_background - summary.n_deg_all;
        wrap.appendChild(el("div", { class: "caption" },
            "Background gene list applied: universe narrowed from " + summary.n_universe_before_background +
            " tested genes to " + summary.n_universe + " (intersected with your uploaded background)." +
            (degsDropped > 0
                ? " " + degsDropped + " DEG(s) outside the background were excluded from every count below."
                : "")));
    }

    if (summary.length_correction_used) {
        wrap.appendChild(alertBox("info", "Gene-length bias correction (Wallenius/GOseq-style) was applied."));
    } else {
        wrap.appendChild(alertBox("warning",
            "Gene-length bias correction was not applied" +
            (summary.length_correction_reason ? ": " + summary.length_correction_reason : ".") +
            " Results use the standard (uncorrected) hypergeometric test."));
    }

    // Per-gene log2FC lookup (uppercase symbol -> log2FC), built once from the
    // retriever's resolved table (not an ORA result table) - feeds the gene
    // volcano and the gene-concept network's gene-node coloring.
    const geneLog2FC = {};
    const resolvedRows = (STATE.transcriptomics.resolved && STATE.transcriptomics.resolved.resolved_table && STATE.transcriptomics.resolved.resolved_table.rows) || [];
    resolvedRows.forEach((r) => {
        const sym = String(r.gene_symbol || "").toUpperCase();
        if (sym && r.log2FC !== null && isFinite(Number(r.log2FC))) geneLog2FC[sym] = Number(r.log2FC);
    });
    const oraCfg = STATE.transcriptomics.oraConfig || {};
    const padjCutoff = typeof oraCfg.padj_cutoff === "number" ? oraCfg.padj_cutoff : 0.05;
    const lfcCutoff = typeof oraCfg.lfc_cutoff === "number" ? oraCfg.lfc_cutoff : 1.0;

    // Significance threshold: display-only (mirrors Genomics' per-database
    // sigThreshold), drives chart filtering, row highlighting, and the
    // hide-charts-entirely gate below. Adjusting it re-renders only the
    // current panel, not the whole screen, so the table selection/scroll is
    // never lost.
    // Assigned further down, once the MSigDB summary exists. The summary sits
    // OUTSIDE panelHost and caches its build, so without this hook it would
    // keep showing bar lengths, colours and an axis title computed at the old
    // alpha while the table right below it used the new one - two
    // contradictory numbers on the same screen.
    let refreshMsigdbSummary = null;

    const alphaRow = el("div", { class: "ora-sig-threshold-row" });
    alphaRow.appendChild(el("span", { class: "chart-option-label" }, "Significance threshold (FDR <)"));
    const alphaInput = el("input", { type: "number", class: "ora-sig-threshold-input", min: "0", max: "1", step: "0.001" });
    alphaInput.value = STATE.transcriptomics.oraSigAlpha;
    alphaInput.addEventListener("change", () => {
        let v = parseFloat(alphaInput.value);
        if (!isFinite(v) || v <= 0 || v > 1) v = 0.05;
        STATE.transcriptomics.oraSigAlpha = v;
        alphaInput.value = v;
        renderPanel(select.value);
        if (refreshMsigdbSummary) refreshMsigdbSummary();
    });
    alphaRow.appendChild(alphaInput);
    wrap.appendChild(alphaRow);

    // Multiple-testing honesty: each database x direction (and, since MSigDB
    // was added, each SELECTED COLLECTION x direction) is independently
    // Benjamini-Hochberg corrected - there is no correction ACROSS these
    // families. Scanning every table and reporting the best hit inflates the
    // true error rate above what any single FDR column states. Same caveat
    // already documented for the Regulon module; sharpened here because
    // selecting several MSigDB collections can push the family count well
    // past what GO+Reactome+one-database ever produced (up to 24).
    if (typeof summary.n_independent_families === "number") {
        wrap.appendChild(el("div", { class: "caption" },
            "This run scanned " + summary.n_independent_families + " independently FDR-corrected " +
            "result table(s) (no correction is applied across them). Selecting more databases/MSigDB " +
            "collections and reporting the best hit across all of them inflates the true error rate " +
            "above what any single FDR column states."));
    }

    // Cross-collection summary: which MSigDB collection is carrying the
    // signal, at a glance, without clicking through every rail entry.
    const hasMsigdbFamilies = Object.keys(tables).some((k) => k.startsWith("MSigDB: ") && k.endsWith(" (all)"));
    if (hasMsigdbFamilies) {
        // Card-wrapped, not a loose button in the page flow - it used to sit
        // directly between two plain-text captions with no border, which
        // read as disconnected from the "Result table" card immediately
        // below it. Matches that card's own header treatment
        // (.result-pane-header/.result-pane-title) so the two read as one
        // consistent stack of sections rather than a button then a card.
        const summaryCard = el("div", { class: "ora-db-card msigdb-summary-card" });
        summaryCard.appendChild(el("div", { class: "result-pane-header" },
            el("h2", { class: "result-pane-title" }, "MSigDB collections")));
        const summaryDiv = el("div", { id: _oraChartId("msigdbsummary", "transcriptomics_ora", null), class: "chart-container", style: { display: "none" } });
        let summaryShown = false, summaryBuilt = false;
        const summaryBtn = el("button", { class: "btn chart-toggle-btn" }, "Show MSigDB collection summary");
        summaryBtn.addEventListener("click", () => {
            if (!summaryShown) {
                summaryDiv.style.display = "";
                if (!summaryBuilt) {
                    buildMsigdbOraCollectionSummary(summaryDiv, tables, STATE.transcriptomics.oraSigAlpha, "transcriptomics_ora");
                    summaryBuilt = true;
                } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(summaryDiv);
                summaryShown = true; summaryBtn.textContent = "Hide MSigDB collection summary";
            } else { summaryDiv.style.display = "none"; summaryShown = false; summaryBtn.textContent = "Show MSigDB collection summary"; }
        });
        // Invalidate on alpha change; rebuild immediately only if it's on
        // screen, otherwise let the next "Show" rebuild it at the new alpha.
        refreshMsigdbSummary = () => {
            summaryBuilt = false;
            if (summaryShown) {
                buildMsigdbOraCollectionSummary(summaryDiv, tables, STATE.transcriptomics.oraSigAlpha, "transcriptomics_ora");
                summaryBuilt = true;
            }
        };
        summaryCard.appendChild(summaryBtn);
        summaryCard.appendChild(summaryDiv);
        wrap.appendChild(summaryCard);
    }

    const sortedKeys = _sortedTranscriptomicsOraKeys(tables);
    const selectorRow = el("div", { class: "row" });
    const select = el("select", { class: "text-input transcriptomics-ora-result-select" });
    sortedKeys.forEach((key) => {
        const block = tables[key];
        const suffix = block.n_terms_tested > 0 ? ` (${block.n_terms_tested} terms)` : " (no terms)";
        select.appendChild(el("option", { value: key }, key + suffix));
    });
    selectorRow.appendChild(el("div", { class: "col" }, select));
    wrap.appendChild(selectorRow);

    const panelHost = el("div", {});
    wrap.appendChild(panelHost);

    function renderPanel(key) {
        panelHost.innerHTML = "";
        const block = tables[key];
        if (!block) return;
        const alpha = STATE.transcriptomics.oraSigAlpha;
        const slug = _transcriptomicsOraSlug(key);
        // Same card chrome as Genomics' renderOraResultBlock (bordered card,
        // colored top border keyed by slug), so the results panel doesn't
        // look visually unfinished next to it.
        const card = el("div", { class: "rank-card ora-db-card ora-result-card ora-db-card-" + slug });
        panelHost.appendChild(card);
        const nav = makeResultNav(card, { stateKey: "transcriptomics_ora" });

        const overviewPane = nav.pane("Overview");
        if (block.note) overviewPane.appendChild(alertBox("info", block.note));
        if (typeof block.n_universe_total === "number" && block.n_terms_tested > 0) {
            overviewPane.appendChild(metricCards([
                ["Terms tested", block.n_terms_tested],
                ["Universe size", block.universe_size],
                ["DEGs in universe", block.n_deg_in_universe],
                ["DEGs not annotated", block.n_deg_total - block.n_deg_in_universe],
            ], 4));
        }
        // "all"/"up"/"down" are refit and corrected independently, so this can
        // legitimately differ from the top summary banner (which reflects
        // "all" only) - e.g. a smaller up/down gene list can yield too flat a
        // length-bias fit even when "all" found a real one.
        if (block.length_correction_used !== summary.length_correction_used) {
            overviewPane.appendChild(el("div", { class: "caption" },
                block.length_correction_used
                    ? "Gene-length bias correction was applied for this specific table."
                    : "Gene-length bias correction was NOT applied for this specific table" +
                      (block.length_correction_reason ? ": " + block.length_correction_reason : ".")));
        }
        if (block.n_terms_excluded_min_size > 0) {
            overviewPane.appendChild(el("div", { class: "caption" },
                block.n_terms_excluded_min_size + " term(s) excluded by the term-size filter."));
        }
        if (block.n_terms_no_overlap > 0) {
            overviewPane.appendChild(el("div", { class: "caption" },
                block.n_terms_no_overlap + " term(s) contained none of the DEGs and were not tested " +
                "(excluded from the FDR calculation, not shown)."));
        }
        if (!block.table.rows || block.table.rows.length === 0) {
            // block.note (already rendered above) already explains WHY the
            // table is empty for every _empty_block()/early-return path (no
            // organism, fetch failure, no DEGs annotated here, etc.) - only
            // fall back to the generic term-size-filter reasons below when
            // there's no such note, otherwise the user sees the correct
            // reason immediately followed by an unrelated, false one.
            const reason = block.note ? null : (
                block.n_terms_no_overlap > 0
                    ? "No term (of the size-eligible ones) contained any of the DEGs for this selection."
                    : "No terms passed the term-size filter for this selection."
            );
            if (reason) overviewPane.appendChild(alertBox("info", reason));
            nav.finish();
            return;
        }

        // Result table - its own pane now, so no separate show/hide toggle is
        // needed: selecting the pane IS the show/hide action.
        const tablePane = nav.pane("Result table");
        tablePane.appendChild(buildTable(block.table, {
            maxHeight: 480,
            rowClassFn: (row) => {
                const fdr = row["FDR"];
                return (fdr !== "" && fdr !== null && isFinite(Number(fdr)) && Number(fdr) < alpha) ? "row-significant" : null;
            },
            cellFormatFn: (col, val) => {
                if (col === "p-value" || col === "FDR" || col === "Bonferroni") return _fmtSci(val, 3);
                if (col === "Fold Enrichment") return _fmtFixed(val, 2);
                return val;
            },
        }));
        const dlBtn = el("button", { class: "btn" }, "Download result table");
        dlBtn.addEventListener("click", () => {
            // Filename must carry the result key: downloadFile() sets
            // a.download, which OVERRIDES the server's own per-key
            // Content-Disposition name. With a fixed name every table in the
            // run saved as the same file (up to 24 of them once MSigDB
            // collections are selected), leaving the user with
            // "..._result.csv", "..._result (1).csv" and no way to tell which
            // family each one was. Slug matches API.py's own re.sub exactly.
            downloadFile("/api/transcriptomics/ora/download?key=" + encodeURIComponent(key),
                         "transcriptomics_ora_" + _transcriptomicsOraSlug(key) + ".csv");
        });
        tablePane.appendChild(dlBtn);

        // Mirrors Genomics' hasSigRows gate: if nothing in THIS table passed
        // the significance threshold, the whole Charts section (including the
        // gene volcano) is hidden rather than showing a degenerate/flat figure.
        const hasSigRows = (block.table.rows || []).some(row =>
            row["FDR"] !== null && row["FDR"] !== "" && isFinite(Number(row["FDR"])) && Number(row["FDR"]) < alpha
        );
        if (!hasSigRows) { nav.finish(); return; }

        function makeChartPane(title) {
            const p = nav.pane(title, { group: "Charts" });
            const controls = el("div", { class: "chart-toggle-row" });
            const area = el("div", { class: "chart-pair" });
            p.appendChild(controls);
            p.appendChild(area);
            const toggles = [];
            return {
                el: p,
                option: (label, options, initialValue, onChange) => {
                    const control = el("label", { class: "chart-option-control" });
                    control.appendChild(el("span", { class: "chart-option-label" }, label));
                    const sel = el("select", { class: "chart-option-select" });
                    options.forEach(([optionLabel, optionValue]) => {
                        sel.appendChild(el("option", { value: optionValue }, optionLabel));
                    });
                    sel.value = initialValue;
                    sel.addEventListener("change", () => onChange(sel.value));
                    control.appendChild(sel);
                    controls.appendChild(control);
                    return sel;
                },
                toggle: (kind, label, buildFn) => {
                    const chartId = _oraChartId(kind, slug, null);
                    const isNetwork = kind === "enrichmap" || kind === "network";
                    const chartDiv = el("div", { id: chartId, class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
                    const chartBtn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
                    let chartShown = false;
                    let chartBuilt = false;
                    function renderChart() { buildFn(chartDiv); chartBuilt = true; }
                    chartBtn.addEventListener("click", () => {
                        if (!chartShown) {
                            chartDiv.style.display = "";
                            if (!chartBuilt) renderChart();
                            else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                            chartShown = true;
                            chartBtn.textContent = "Hide " + label;
                        } else {
                            chartDiv.style.display = "none";
                            chartShown = false;
                            chartBtn.textContent = "Show " + label;
                        }
                    });
                    controls.appendChild(chartBtn);
                    area.appendChild(chartDiv);
                    const t = { redraw: () => { if (chartShown) renderChart(); } };
                    toggles.push(t);
                    return t;
                },
                redraw: () => toggles.forEach(t => t.redraw()),
            };
        }

        // --- Gene set overview: dot, bar, gene volcano ----------------------
        const overview = makeChartPane("Gene set overview");
        overview.el.insertBefore(el("div", { class: "ora-chart-threshold-note" },
            "Significance threshold: FDR < " + alpha + " (used for table highlighting, chart filtering, and gene-volcano coloring; adjustable above)."),
            overview.el.firstChild);
        let dotColorScale1 = ORA_DOT_COLOR_SCALES[0][1];
        let barColor = ORA_BAR_COLORS[0][1];
        let topN = 20;
        let sortBy = "fdr";
        let fontSize = 12;

        // Terms shown/Rank & show by are shared with the relationships group
        // below (enrichment map + gene-concept network also plot by topN/
        // sortBy) - both panes must redraw, not just this one.
        overview.option("Terms shown", [
            ["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"],
            ["50", "50"], ["100", "100"], ["200", "200"], ["All significant", "all"],
        ], String(topN === Infinity ? "all" : topN), (value) => {
            topN = value === "all" ? Infinity : Number(value);
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Rank & show by", [["FDR", "fdr"], ["Fold enrichment", "fold"], ["Gene count", "k"], ["Gene ratio", "ratio"]], sortBy, (value) => {
            sortBy = value;
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], String(fontSize), (value) => {
            fontSize = Number(value);
            overview.redraw();
            relationships.redraw();
        });

        const dotToggle = overview.toggle("dot", "dot plot", (div) => buildDotPlot(div, block.table, slug, null, dotColorScale1, { topN, sortBy, fontSize, alpha }));
        overview.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale1, (value) => {
            dotColorScale1 = value;
            dotToggle.redraw();
        });
        const barToggle = overview.toggle("bar", "bar plot", (div) => buildBarPlot(div, block.table, slug, null, barColor, { topN, sortBy, fontSize, alpha }));
        overview.option("Bar color", ORA_BAR_COLORS, barColor, (value) => {
            barColor = value;
            barToggle.redraw();
        });
        // Gene-level volcano (transcriptomics-specific): reads the retriever's
        // resolved table directly, not block.table, so it is unaffected by
        // Terms shown/Rank & show by/color controls above (those are term-
        // enrichment concepts; the volcano is a per-gene view of the whole
        // experiment).
        overview.toggle("genevolcano", "gene volcano", (div) =>
            buildTranscriptomicsGeneVolcano(div, resolvedRows, slug, padjCutoff, lfcCutoff, { fontSize }));

        // --- Gene set relationships: enrichment map, gene-concept network ---
        const relationships = makeChartPane("Gene set relationships");
        relationships.el.insertBefore(el("div", { class: "ora-chart-threshold-note ora-chart-encoding-note" },
            "Circle color = -log10(FDR) (darker = more significant); circle size = number of DEGs in the term. " +
            "Rank & show by and Terms shown (Gene set overview) choose which terms are drawn."),
            relationships.el.firstChild);
        let dotColorScale2 = ORA_DOT_COLOR_SCALES[0][1];
        let edgeThreshold = 1;
        const enrichMapToggle = relationships.toggle("enrichmap", "enrichment map", (div) =>
            buildEnrichmentMap(div, block.table, slug, null, dotColorScale2, { topN, sortBy, fontSize, edgeThreshold, alpha }));
        relationships.option("Edge threshold (shared genes)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (value) => {
            edgeThreshold = Number(value);
            enrichMapToggle.redraw();
        });
        const networkToggle = relationships.toggle("network", "gene-concept network", (div) =>
            buildTranscriptomicsGeneConceptNetwork(div, block.table, slug, geneLog2FC, { topN, sortBy, fontSize, alpha, dotColorScale: dotColorScale2 }));
        relationships.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale2, (value) => {
            dotColorScale2 = value;
            enrichMapToggle.redraw();
            networkToggle.redraw();
        });

        nav.finish();
    }

    select.addEventListener("change", () => renderPanel(select.value));
    if (sortedKeys.length > 0) {
        select.value = sortedKeys[0];
        renderPanel(sortedKeys[0]);
    }

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.transcriptomics.oraConfig = null;
        STATE.transcriptomics.oraResults = null;
        STATE.transcriptomics.oraRunLocked = false;
        STATE.transcriptomics._oraDraft = null;
        STATE.aaPage = "transcriptomics_ora_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// Proteomics — tile grid, upload/resolve, Over-Representation Analysis
//
// Identifier space is UniProt accession ONLY - no identifierType choice,
// unlike Transcriptomics (see Proteomics_module_retriever.py's module
// docstring for why). Only ORA is wired this milestone; the other three
// tiles (Quantitative Enrichment, PPI Network, PTM) are visible but
// disabled, matching TRANSCRIPTOMICS_TILES_NOT_YET_IMPLEMENTED's pattern.
// ---------------------------------------------------------------------------

const PROTEOMICS_TILES_NOT_YET_IMPLEMENTED = new Set([
    "Quantitative Enrichment", "PPI Network / Module Detection", "PTM / Phosphosite Enrichment",
]);

// Mirrors Proteomics_module.PROTEOMICS_ORA_DATABASES - keep in sync.
const PROTEOMICS_ORA_DATABASES = [
    "GO", "Reactome", "MSigDB", "UniProt (disease)", "UniProt (keywords)",
    "InterPro", "Complex Portal", "DISEASES", "Custom database",
];
// DISEASES (Jensen Lab) is a HUMAN-only resource by design, so its
// checkbox greys out for any other organism - same treatment MSigDB gets
// for being human/mouse-only, rather than letting a run come back
// silently empty with no explanation.
function _diseasesAvailableForTaxon(taxonId) {
    return String(taxonId) === "9606";
}

function renderProteomicsScreen() {
    const wrap = el("div", {});
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle" }, "Proteomics"));
    const grid = el("div", { class: "aa-genomics-grid" });
    [
        "Over-Representation Analysis",
        "Quantitative Enrichment",
        "PPI Network / Module Detection",
        "PTM / Phosphosite Enrichment",
    ].forEach((name) => {
        const notYetImplemented = PROTEOMICS_TILES_NOT_YET_IMPLEMENTED.has(name);
        const btn = el("button", {
            class: "module-btn proteomics" + (notYetImplemented ? " module-btn-disabled" : ""),
            title: notYetImplemented ? name + " is not implemented yet - coming soon." : "",
        }, name + (notYetImplemented ? " (coming soon)" : ""));
        btn.disabled = notYetImplemented;
        btn.addEventListener("click", () => {
            if (notYetImplemented) return;
            STATE.aaPage = "proteomics_module";
            STATE.proteomicsModuleName = name;
            STATE.proteomics = freshProteomicsState();
            render();
        });
        grid.appendChild(btn);
    });
    wrap.appendChild(grid);
    return wrap;
}

function backToProteomicsActions() {
    const actions = el("div", { class: "right-actions" });
    const btn = el("button", { class: "btn" }, "Back to Proteomics");
    btn.addEventListener("click", () => { STATE.aaPage = "proteomics"; render(); });
    actions.appendChild(el("div", { class: "btn-row", style: { width: "25%", minWidth: "220px" } }, btn));
    return actions;
}

function renderProteomicsOrganismSelector() {
    const wrap = el("div", {});
    const input = el("input", { class: "text-input", type: "text",
        placeholder: "Species / organism (name or UniProt taxon ID)" });
    const resultsSelect = el("select", { class: "hidden" });
    const selectedCaption = el("div", { class: "caption" });

    function refreshCaption() {
        selectedCaption.textContent = STATE.proteomics.taxonId
            ? "Selected: " + STATE.proteomics.taxonName + " (taxon " + STATE.proteomics.taxonId + ")"
            : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = await apiJson("/api/contexts/search?database=" + encodeURIComponent("Rhea + UniProt") + "&query=" + encodeURIComponent(query) + "&limit=10");
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                (data.results || []).forEach((rec, idx) => {
                    const label = rec.context_id && rec.context_name ? rec.context_id + " | " + rec.context_name : (rec.context_name || rec.context_id || "");
                    resultsSelect.appendChild(el("option", { value: String(idx) }, label));
                });
                resultsSelect._records = data.results || [];
                resultsSelect.classList.toggle("hidden", (data.results || []).length === 0);
            } catch (err) { /* keep prior dropdown contents on transient search errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        const newTaxonId = String(rec.context_id || "");
        const hadResolvedForDifferentOrganism =
            !!STATE.proteomics.resolved && newTaxonId !== STATE.proteomics.taxonId;
        STATE.proteomics.taxonId = newTaxonId;
        STATE.proteomics.taxonName = String(rec.context_name || "");
        if (hadResolvedForDifferentOrganism) {
            STATE.proteomics.resolved = null;
            STATE.proteomics.oraConfig = null;
            STATE.proteomics.oraResults = null;
            STATE.proteomics.oraRunLocked = false;
        }
        render();
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function renderProteomicsUploadSection() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "section-label" }, "Upload protein-level differential-abundance table"));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.proteomics.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            try {
                const data = await apiUpload("/api/proteomics/upload", formData);
                STATE.proteomics.deTable = data.preview;
                STATE.proteomics.uploadedFilename = file.name;
                STATE.proteomics.columnMapping = data.column_mapping;
                // Server's upload endpoint resets the whole session sub-dict -
                // clear client-side too, or the UI would keep showing state
                // the server can no longer serve (mirrors Transcriptomics).
                STATE.proteomics.customDb = null;
                STATE.proteomics.background = null;
                STATE.proteomics.resolved = null;
                STATE.proteomics.oraConfig = null;
                STATE.proteomics.oraResults = null;
                STATE.proteomics.oraRunLocked = false;
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));

    wrap.appendChild(el("div", { class: "caption" },
        "Required columns: UniProt accession (or protein group, e.g. \"P04637;P04637-2\") + log2FC + " +
        "pvalue. Optional: padj/FDR (computed via Benjamini-Hochberg if omitted); MaxQuant's QC flag " +
        "columns (\"Potential contaminant\", \"Reverse\", \"Only identified by site\") are auto-detected " +
        "and those rows removed, alongside any CON__/REV__-prefixed ID."));

    const exLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Example file");
    exLink.addEventListener("click", () => downloadFile("/api/proteomics/example", "example_proteomics_de_full.csv"));
    wrap.appendChild(exLink);
    const messyLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Messy MaxQuant-style example");
    messyLink.addEventListener("click", () => downloadFile("/api/proteomics/example?messy=true", "example_proteomics_maxquant_messy.csv"));
    wrap.appendChild(document.createTextNode(" "));
    wrap.appendChild(messyLink);
    wrap.appendChild(el("div", { class: "caption" },
        "The messy example deliberately contains protein groups, CON__/REV__ rows, and QC flag columns - " +
        "upload it to see the contaminant-filtering report below."));

    if (STATE.proteomics.deTable) {
        if (STATE.proteomics.columnMapping) {
            const cm = STATE.proteomics.columnMapping;
            const parts = Object.keys(cm).map((k) => k + " -> " + cm[k]);
            wrap.appendChild(el("div", { class: "caption" }, "Detected columns: " + parts.join(", ")));
        }
        // Dynamic height matching renderTranscriptomicsUploadSection's own
        // sizing exactly - a fixed 300px was an unjustified inconsistency
        // (also: the preview itself now carries up to 200 rows, not 10).
        const maxHeight = Math.min(400, Math.max(200, (STATE.proteomics.deTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.proteomics.deTable, { maxHeight }));
    }

    return wrap;
}

function renderProteomicsResolveControls() {
    const wrap = el("div", {});
    const disabled = !STATE.proteomics.deTable || !STATE.proteomics.taxonId;
    const btn = el("button", { class: "btn" }, "Resolve accessions");
    btn.disabled = disabled;
    btn.addEventListener("click", async () => {
        btn.disabled = true;
        const spinner = el("div", { class: "caption" }, "Resolving accessions via UniProt ID Mapping...");
        wrap.appendChild(spinner);
        try {
            const { job_id } = await apiPostJson("/api/proteomics/resolve", {
                taxon_id: STATE.proteomics.taxonId,
                taxon_name: STATE.proteomics.taxonName,
            });
            pollJob(job_id, {
                onDone: (result) => {
                    STATE.proteomics.resolved = result;
                    render();
                },
                onError: (err) => { alert(err); render(); },
            });
        } catch (err) {
            alert(err.message);
            render();
        }
    });
    wrap.appendChild(btn);
    return wrap;
}

function renderProteomicsResults() {
    const wrap = el("div", {});
    const resolved = STATE.proteomics.resolved;
    const meta = resolved.metrics;
    wrap.appendChild(el("div", { class: "section-label" }, "Resolved protein table"));
    wrap.appendChild(buildTable(resolved.resolved_table, { maxHeight: 420 }));
    wrap.appendChild(metricCards([
        ["Mapped", meta.n_mapped],
        ["Unmapped", meta.n_unmapped],
        ["Resolved via secondary accession", meta.n_resolved_secondary],
        ["Duplicates collapsed", meta.n_duplicates_collapsed],
        ["Dropped (non-numeric)", meta.n_dropped_nonnumeric],
    ], 5));

    // Contaminant/decoy report - the proteomics-specific behaviour with no
    // transcriptomics equivalent, so it earns its own visible block rather
    // than being buried in the generic metric cards above.
    const c = meta.contaminants || {};
    if (c.n_total_removed > 0) {
        wrap.appendChild(alertBox("info",
            c.n_total_removed + " contaminant/decoy row(s) removed before resolution: " +
            c.n_prefix_removed + " by CON__/REV__ ID prefix" +
            (c.flag_columns_found && c.flag_columns_found.length
                ? ", " + c.n_flagged_removed + " by QC flag column(s) (" + c.flag_columns_found.join(", ") + ")"
                : "") + "."));
    }
    if (meta.padj_computed) {
        wrap.appendChild(alertBox("info",
            "No adjusted p-value (padj/FDR) column was found - it was computed automatically via " +
            "Benjamini-Hochberg from the p-value column."));
    }
    if (meta.n_fetch_failures > 0) {
        wrap.appendChild(alertBox("warning",
            meta.n_fetch_failures + " UniProt ID Mapping batch(es) failed after retries. Consider " +
            "re-running if this table looks incomplete."));
    }

    const dlResolvedBtn = el("button", { class: "btn" }, "Download resolved table (CSV)");
    dlResolvedBtn.addEventListener("click", () => {
        downloadFile("/api/proteomics/download/resolved.csv", "proteomics_resolved.csv");
    });
    wrap.appendChild(dlResolvedBtn);

    if (STATE.proteomicsModuleName === "Over-Representation Analysis") {
        const continueBtn = el("button", { class: "btn btn-primary", style: { marginTop: "0.8rem" } }, "Continue to Over-Representation Analysis");
        continueBtn.addEventListener("click", () => {
            STATE.proteomics.oraConfig = null;
            STATE.proteomics.oraResults = null;
            STATE.proteomics.oraRunLocked = false;
            STATE.aaPage = "proteomics_ora_config";
            render();
        });
        wrap.appendChild(continueBtn);
    } else {
        wrap.appendChild(el("div", { class: "caption" }, "This analysis module is not available yet."));
    }
    return wrap;
}

function renderProteomicsModuleScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Proteomics — " + STATE.proteomicsModuleName));
    wrap.appendChild(el("div", { class: "section-label" }, "Organism"));
    wrap.appendChild(renderProteomicsOrganismSelector());
    wrap.appendChild(renderProteomicsUploadSection());
    wrap.appendChild(renderProteomicsResolveControls());
    if (STATE.proteomics.resolved) wrap.appendChild(renderProteomicsResults());
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToProteomicsActions());
    return wrap;
}

function renderProteomicsCustomDbWidget(locked) {
    const wrap = el("div", { class: "ora-custom-db-widget" });
    wrap.appendChild(el("div", { class: "section-label" }, "Custom database"));
    wrap.appendChild(el("div", { class: "caption" },
        "Upload your own gene sets (GMT or a tidy term/gene table), in UniProt-accession space."));

    if (!locked) {
        wrap.appendChild(fileInputWidget({
            accept: ".gmt,.txt,.tsv,.csv",
            buttonLabel: "Choose gene-set file",
            initialFilename: STATE.proteomics.customDb ? STATE.proteomics.customDb.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                try {
                    const data = await apiUpload("/api/proteomics/ora/custom-db-upload", formData);
                    STATE.proteomics.customDb = data;
                    render();
                } catch (err) {
                    STATE.proteomics.customDb = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    if (STATE.proteomics.customDb) {
        const s = STATE.proteomics.customDb;
        wrap.appendChild(el("div", { class: "caption ora-custom-summary" },
            "Loaded \"" + s.filename + "\": " + s.n_terms + " gene set(s), " + s.universe_size + " unique accessions."));
    }

    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "(example input files)");
    exampleLink.addEventListener("click", () => {
        downloadFile("/api/proteomics/ora/custom-db-example", "proteomics_custom_geneset_db_example.zip");
    });
    wrap.appendChild(exampleLink);
    return wrap;
}

function renderProteomicsBackgroundWidget(locked, draft) {
    const wrap = el("div", { class: "ora-background-widget" });

    const toggleRow = el("div", { class: "ora-iea-row ora-background-toggle-row" });
    const toggleCb = el("input", { type: "checkbox", class: "ora-iea-checkbox ora-background-checkbox" });
    toggleCb.checked = draft.useBackground;
    toggleCb.disabled = locked;
    toggleCb.addEventListener("change", () => { draft.useBackground = toggleCb.checked; render(); });
    toggleRow.appendChild(toggleCb);
    toggleRow.appendChild(document.createTextNode(" Background accession list (optional)"));
    wrap.appendChild(toggleRow);

    if (!draft.useBackground) return wrap;

    const body = el("div", { class: "ora-background-body" });
    if (!locked) {
        body.appendChild(fileInputWidget({
            accept: ".txt,.csv,.tsv",
            buttonLabel: "Choose background file",
            initialFilename: STATE.proteomics.background ? STATE.proteomics.background.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                try {
                    const data = await apiUpload("/api/proteomics/ora/background-upload", formData);
                    STATE.proteomics.background = { filename: data.filename, n_accessions: data.n_accessions };
                    render();
                } catch (err) {
                    STATE.proteomics.background = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    if (STATE.proteomics.background) {
        const s = STATE.proteomics.background;
        body.appendChild(el("div", { class: "caption ora-background-summary" },
            "Background: \"" + s.filename + "\" - " + s.n_accessions + " UniProt accession(s). Restricts " +
            "the universe (for every selected database) to accessions present in both this list and the " +
            "proteins quantified in this experiment."));
        if (!locked) {
            const removeBtn = el("button", { class: "btn btn-small ora-background-remove", type: "button" }, "Remove background");
            removeBtn.addEventListener("click", async () => {
                try { await apiPostJson("/api/proteomics/ora/background-remove", {}); } catch (err) { /* ignore */ }
                STATE.proteomics.background = null;
                render();
            });
            body.appendChild(removeBtn);
        }
    } else {
        const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "example file");
        exampleLink.addEventListener("click", () => {
            downloadFile("/api/proteomics/ora/background-example", "example_proteomics_background.csv");
        });
        body.appendChild(el("div", { class: "caption ora-background-help" }, [
            "Provide UniProt accessions (one per line, or a column headed \"UniProt Accession\") for a " +
            "narrower reference set. Leave off to use the default: every protein quantified in this " +
            "experiment. (",
            exampleLink,
            ")",
        ]));
    }
    wrap.appendChild(body);
    return wrap;
}

function renderProteomicsOraConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Proteomics — Over-Representation Analysis: database"));

    const locked = STATE.proteomics.oraRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    if (!STATE.proteomics._oraDraft) {
        STATE.proteomics._oraDraft = {
            databases: [], msigdbCollections: [],
            padjCutoff: 0.05, lfcCutoff: 1.0,
            minTermSize: 10, maxTermSize: 500,
            includeIea: true, useBackground: false,
        };
    }
    const draft = STATE.proteomics._oraDraft;

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const runBtnWrap = el("div", { class: "col" });

    function renderOraRunBtn() {
        runBtnWrap.innerHTML = "";
        const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
        const needsCustomDb = draft.databases.indexOf("Custom database") !== -1 && !STATE.proteomics.customDb;
        const needsMsigdbCollection = draft.databases.indexOf("MSigDB") !== -1 && draft.msigdbCollections.length === 0;
        runBtn.disabled = locked || draft.databases.length === 0 || needsCustomDb || needsMsigdbCollection;
        runBtn.addEventListener("click", () => startProteomicsOraRun(draft));
        runBtnWrap.appendChild(runBtn);
    }

    const taxonId = STATE.proteomics.taxonId;
    const msigdbAvailable = _msigdbAvailableForTaxon(taxonId);
    const diseasesAvailable = _diseasesAvailableForTaxon(taxonId);
    PROTEOMICS_ORA_DATABASES.forEach((dbName) => {
        const unavailable = (dbName === "MSigDB" && !msigdbAvailable)
                          || (dbName === "DISEASES" && !diseasesAvailable);
        if (unavailable) draft.databases = draft.databases.filter((d) => d !== dbName);
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked || unavailable;
        cb.checked = !unavailable && draft.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.databases.indexOf(dbName) === -1) draft.databases.push(dbName);
            } else {
                draft.databases = draft.databases.filter((d) => d !== dbName);
            }
            if (dbName === "Custom database" || dbName === "MSigDB") render();
            else renderOraRunBtn();
        });
        label.appendChild(cb);
        const unavailableNote = dbName === "DISEASES" ? " (human only)" : " (human/mouse only)";
        label.appendChild(document.createTextNode(" " + dbName + (unavailable ? unavailableNote : "")));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    if (!msigdbAvailable) {
        wrap.appendChild(el("div", { class: "caption" },
            "MSigDB is disabled: it covers human and mouse only."));
    }
    if (!diseasesAvailable) {
        wrap.appendChild(el("div", { class: "caption" },
            "DISEASES is disabled: it is a human-only disease-gene resource."));
    }
    if (draft.databases.indexOf("MSigDB") !== -1 && msigdbAvailable) {
        wrap.appendChild(renderMsigdbCollectionSelector(locked, draft, taxonId, renderOraRunBtn));
    }
    if (draft.databases.indexOf("Custom database") !== -1) {
        wrap.appendChild(renderProteomicsCustomDbWidget(locked));
    }
    wrap.appendChild(el("div", { class: "caption" },
        "GO / Reactome / MSigDB / Custom: same engines as Transcriptomics. UniProt (disease) and " +
        "UniProt (keywords): UniProt's own curated annotations - note UniProt's disease entries are " +
        "overwhelmingly one-gene rare diseases, so most are too small to test. DISEASES (Jensen Lab, " +
        "CC BY 4.0, human only): curated gene-disease associations grouped into broad disease " +
        "CATEGORIES (e.g. \"Neurodegenerative disease\", \"Cancer\"), which is what makes it usable for " +
        "enrichment where UniProt's disease entries are not. InterPro: protein domain/family " +
        "membership. Complex Portal: known macromolecular complexes (CC0). Each is tested against all " +
        "significant proteins, and separately against up- and down-regulated ones - independently " +
        "FDR-corrected, never pooled."));

    wrap.appendChild(renderProteomicsBackgroundWidget(locked, draft));

    function numberField(labelText, value, step, onChange) {
        const col = el("div", { class: "col" });
        col.appendChild(el("div", { class: "chart-option-label" }, labelText));
        const input = el("input", { type: "number", class: "text-input", step: step });
        input.value = value;
        input.disabled = locked;
        input.addEventListener("change", () => {
            const v = parseFloat(input.value);
            if (isFinite(v)) onChange(v);
        });
        col.appendChild(input);
        return col;
    }

    const thresholdRow = el("div", { class: "row", style: { marginTop: "0.8rem" } });
    thresholdRow.appendChild(numberField("Adjusted p-value (padj) cutoff", draft.padjCutoff, "0.001", (v) => { draft.padjCutoff = v; }));
    thresholdRow.appendChild(numberField("|log2 Fold-Change| cutoff", draft.lfcCutoff, "0.1", (v) => { draft.lfcCutoff = v; }));
    wrap.appendChild(thresholdRow);

    const termSizeRow = el("div", { class: "row" });
    const minCol = el("div", { class: "col" });
    minCol.appendChild(el("div", { class: "chart-option-label" }, "Minimum term size (proteins)"));
    const minInput = el("input", { type: "number", class: "text-input", min: "1", step: "1" });
    minInput.value = draft.minTermSize;
    minInput.disabled = locked;
    minInput.addEventListener("change", () => {
        const v = parseInt(minInput.value, 10);
        draft.minTermSize = (isFinite(v) && v >= 1) ? v : 10;
        minInput.value = draft.minTermSize;
    });
    minCol.appendChild(minInput);
    termSizeRow.appendChild(minCol);

    const maxCol = el("div", { class: "col" });
    maxCol.appendChild(el("div", { class: "chart-option-label" }, "Maximum term size (proteins, blank = no limit)"));
    const maxInput = el("input", { type: "number", class: "text-input", min: "1", step: "1", placeholder: "no limit" });
    maxInput.value = draft.maxTermSize;
    maxInput.disabled = locked;
    maxInput.addEventListener("change", () => {
        const raw = maxInput.value.trim();
        draft.maxTermSize = raw === "" ? "" : parseInt(raw, 10);
    });
    maxCol.appendChild(maxInput);
    termSizeRow.appendChild(maxCol);
    wrap.appendChild(termSizeRow);

    const ieaRow = el("div", { class: "row", style: { alignItems: "center", marginTop: "0.5rem" } });
    const ieaLabel = el("label", { style: { display: "flex", alignItems: "center", gap: "0.4rem" } });
    const ieaCheckbox = el("input", { type: "checkbox" });
    ieaCheckbox.checked = draft.includeIea;
    ieaCheckbox.disabled = locked;
    ieaCheckbox.addEventListener("change", () => { draft.includeIea = ieaCheckbox.checked; });
    ieaLabel.appendChild(ieaCheckbox);
    ieaLabel.appendChild(document.createTextNode(" Include electronic annotations (IEA)"));
    ieaRow.appendChild(ieaLabel);
    wrap.appendChild(ieaRow);
    wrap.appendChild(el("div", { class: "caption" },
        "Applies to GO and Reactome. Uncheck to restrict both to curated (non-electronic) evidence only."));

    const actions = el("div", { class: "right-actions", style: { marginTop: "1rem" } });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "proteomics_module"; render(); });
    btnRow.appendChild(backBtn);
    renderOraRunBtn();
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function startProteomicsOraRun(draft) {
    if (draft.databases.length === 0) return;
    const payload = {
        msigdb_collections: draft.databases.indexOf("MSigDB") !== -1 ? draft.msigdbCollections : [],
        databases: draft.databases,
        padj_cutoff: draft.padjCutoff,
        lfc_cutoff: draft.lfcCutoff,
        min_term_size: draft.minTermSize,
        max_term_size: draft.maxTermSize === "" ? null : draft.maxTermSize,
        include_iea: draft.includeIea,
        use_background: !!(draft.useBackground && STATE.proteomics.background),
    };
    apiPostJson("/api/proteomics/ora/configure", payload).then(() => {
        STATE.proteomics.oraConfig = payload;
        STATE.proteomics.oraRunLocked = true;
        STATE.proteomics.oraResults = null;
        STATE.aaPage = "proteomics_ora_results";
        render();
        return apiPostJson("/api/proteomics/ora/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("proteomics-ora-progress-inner");
                const label = document.getElementById("proteomics-ora-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => { STATE.proteomics.oraResults = result; render(); },
            onError: (err) => {
                alert("ORA run failed: " + err);
                STATE.proteomics.oraRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.proteomics.oraRunLocked = false;
        render();
    });
}

const PROTEOMICS_ORA_FAMILY_ORDER = [
    "GO: Biological Process", "GO: Molecular Function", "GO: Cellular Component",
    "Reactome", "UniProt disease", "UniProt keywords", "InterPro", "Complex Portal",
    "DISEASES",
    ...TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.map((c) => "MSigDB: " + c.label),
    "Custom database",
];
const PROTEOMICS_ORA_DIRECTION_ORDER = ["all", "up", "down"];

function _splitProteomicsOraKey(key) {
    const m = key.match(/^(.*) \((all|up|down)\)$/);
    return m ? [m[1], m[2]] : [key, "all"];
}
function _sortedProteomicsOraKeys(tables) {
    return Object.keys(tables).sort((a, b) => {
        const [famA, dirA] = _splitProteomicsOraKey(a);
        const [famB, dirB] = _splitProteomicsOraKey(b);
        const iA = PROTEOMICS_ORA_FAMILY_ORDER.indexOf(famA);
        const iB = PROTEOMICS_ORA_FAMILY_ORDER.indexOf(famB);
        if (iA !== iB) return iA - iB;
        return PROTEOMICS_ORA_DIRECTION_ORDER.indexOf(dirA) - PROTEOMICS_ORA_DIRECTION_ORDER.indexOf(dirB);
    });
}
function _proteomicsOraSlug(key) {
    return key.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "");
}

function renderProteomicsOraResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Proteomics — Over-Representation Analysis: results"));

    if (!STATE.proteomics.oraResults) {
        wrap.appendChild(alertBox("info", "Running analysis..."));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "proteomics-ora-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "proteomics-ora-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const { tables, summary } = STATE.proteomics.oraResults;

    wrap.appendChild(metricCards([
        ["Universe (quantified proteins)", summary.n_universe],
        ["Significant (all)", summary.n_sig_all],
        ["Significant (up)", summary.n_sig_up],
        ["Significant (down)", summary.n_sig_down],
    ], 4));

    if (summary.background_applied) {
        const sigDropped = summary.n_sig_all_before_background - summary.n_sig_all;
        wrap.appendChild(el("div", { class: "caption" },
            "Background accession list applied: universe narrowed from " + summary.n_universe_before_background +
            " quantified proteins to " + summary.n_universe + " (intersected with your uploaded background)." +
            (sigDropped > 0
                ? " " + sigDropped + " significant protein(s) outside the background were excluded from every count below."
                : "")));
    }

    // Per-protein log2FC lookup (accession -> log2FC), built once from the
    // retriever's resolved table - feeds the protein volcano and the
    // protein-concept network's node coloring (reuses Transcriptomics'
    // builders unchanged: both key off table cell TEXT, which is accessions
    // here rather than gene symbols, and the label fallback in the volcano
    // builder already reads gene_symbol/input_id from these exact rows).
    const proteinLog2FC = {};
    const resolvedRows = (STATE.proteomics.resolved && STATE.proteomics.resolved.resolved_table && STATE.proteomics.resolved.resolved_table.rows) || [];
    resolvedRows.forEach((r) => {
        const acc = String(r.accession || "").toUpperCase();
        if (acc && r.log2FC !== null && isFinite(Number(r.log2FC))) proteinLog2FC[acc] = Number(r.log2FC);
    });
    const oraCfg = STATE.proteomics.oraConfig || {};
    const padjCutoff = typeof oraCfg.padj_cutoff === "number" ? oraCfg.padj_cutoff : 0.05;
    const lfcCutoff = typeof oraCfg.lfc_cutoff === "number" ? oraCfg.lfc_cutoff : 1.0;

    let refreshMsigdbSummary = null;

    const alphaRow = el("div", { class: "ora-sig-threshold-row" });
    alphaRow.appendChild(el("span", { class: "chart-option-label" }, "Significance threshold (FDR <)"));
    const alphaInput = el("input", { type: "number", class: "ora-sig-threshold-input", min: "0", max: "1", step: "0.001" });
    alphaInput.value = STATE.proteomics.oraSigAlpha;
    alphaInput.addEventListener("change", () => {
        let v = parseFloat(alphaInput.value);
        if (!isFinite(v) || v <= 0 || v > 1) v = 0.05;
        STATE.proteomics.oraSigAlpha = v;
        alphaInput.value = v;
        renderPanel(select.value);
        if (refreshMsigdbSummary) refreshMsigdbSummary();
    });
    alphaRow.appendChild(alphaInput);
    wrap.appendChild(alphaRow);

    wrap.appendChild(el("div", { class: "caption" },
        "Each database x direction (and each selected MSigDB collection x direction) is independently " +
        "Benjamini-Hochberg corrected - there is no correction ACROSS these families. Scanning every " +
        "table and reporting the best hit across all of them inflates the true error rate above what " +
        "any single FDR column states."));

    const hasMsigdbFamilies = Object.keys(tables).some((k) => k.startsWith("MSigDB: ") && k.endsWith(" (all)"));
    if (hasMsigdbFamilies) {
        const summaryCard = el("div", { class: "ora-db-card msigdb-summary-card" });
        summaryCard.appendChild(el("div", { class: "result-pane-header" },
            el("h2", { class: "result-pane-title" }, "MSigDB collections")));
        const summaryDiv = el("div", { id: _oraChartId("msigdbsummary", "proteomics_ora", null), class: "chart-container", style: { display: "none" } });
        let summaryShown = false, summaryBuilt = false;
        const summaryBtn = el("button", { class: "btn chart-toggle-btn" }, "Show MSigDB collection summary");
        summaryBtn.addEventListener("click", () => {
            if (!summaryShown) {
                summaryDiv.style.display = "";
                if (!summaryBuilt) {
                    buildMsigdbOraCollectionSummary(summaryDiv, tables, STATE.proteomics.oraSigAlpha, "proteomics_ora");
                    summaryBuilt = true;
                } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(summaryDiv);
                summaryShown = true; summaryBtn.textContent = "Hide MSigDB collection summary";
            } else { summaryDiv.style.display = "none"; summaryShown = false; summaryBtn.textContent = "Show MSigDB collection summary"; }
        });
        refreshMsigdbSummary = () => {
            summaryBuilt = false;
            if (summaryShown) {
                buildMsigdbOraCollectionSummary(summaryDiv, tables, STATE.proteomics.oraSigAlpha, "proteomics_ora");
                summaryBuilt = true;
            }
        };
        summaryCard.appendChild(summaryBtn);
        summaryCard.appendChild(summaryDiv);
        wrap.appendChild(summaryCard);
    }

    const sortedKeys = _sortedProteomicsOraKeys(tables);
    const selectorRow = el("div", { class: "row" });
    const select = el("select", { class: "text-input transcriptomics-ora-result-select" });
    sortedKeys.forEach((key) => {
        const block = tables[key];
        const suffix = block.n_terms_tested > 0 ? ` (${block.n_terms_tested} terms)` : " (no terms)";
        select.appendChild(el("option", { value: key }, key + suffix));
    });
    selectorRow.appendChild(el("div", { class: "col" }, select));
    wrap.appendChild(selectorRow);

    const panelHost = el("div", {});
    wrap.appendChild(panelHost);

    function renderPanel(key) {
        panelHost.innerHTML = "";
        const block = tables[key];
        if (!block) return;
        const alpha = STATE.proteomics.oraSigAlpha;
        const slug = _proteomicsOraSlug(key);
        const card = el("div", { class: "rank-card ora-db-card ora-result-card ora-db-card-" + slug });
        panelHost.appendChild(card);
        const nav = makeResultNav(card, { stateKey: "proteomics_ora" });

        const overviewPane = nav.pane("Overview");
        if (block.caveat) overviewPane.appendChild(alertBox("warning", block.caveat));
        if (block.note) overviewPane.appendChild(alertBox("info", block.note));
        if (block.n_terms_tested > 0) {
            overviewPane.appendChild(metricCards([
                ["Terms tested", block.n_terms_tested],
                ["Universe size", block.universe_size],
                ["Significant proteins in universe", block.n_input_in_universe],
                ["Excluded by term-size filter", block.n_terms_excluded_min_size],
            ], 4));
        }
        if (block.n_terms_no_overlap > 0) {
            overviewPane.appendChild(el("div", { class: "caption" },
                block.n_terms_no_overlap + " term(s) contained none of the significant proteins and were " +
                "not tested (excluded from the FDR calculation, not shown)."));
        }
        if (!block.table.rows || block.table.rows.length === 0) {
            const reason = block.note || block.caveat ? null : (
                block.n_terms_no_overlap > 0
                    ? "No term (of the size-eligible ones) contained any significant protein for this selection."
                    : "No terms passed the term-size filter for this selection."
            );
            if (reason) overviewPane.appendChild(alertBox("info", reason));
            nav.finish();
            return;
        }

        const tablePane = nav.pane("Result table");
        tablePane.appendChild(buildTable(block.table, {
            maxHeight: 480,
            rowClassFn: (row) => {
                const fdr = row["FDR"];
                return (fdr !== "" && fdr !== null && isFinite(Number(fdr)) && Number(fdr) < alpha) ? "row-significant" : null;
            },
            cellFormatFn: (col, val) => {
                if (col === "p-value" || col === "FDR" || col === "Bonferroni") return _fmtSci(val, 3);
                if (col === "Fold Enrichment") return _fmtFixed(val, 2);
                return val;
            },
        }));
        const dlBtn = el("button", { class: "btn" }, "Download result table");
        dlBtn.addEventListener("click", () => {
            downloadFile("/api/proteomics/ora/download?key=" + encodeURIComponent(key),
                         "proteomics_ora_" + _proteomicsOraSlug(key) + ".csv");
        });
        tablePane.appendChild(dlBtn);

        const hasSigRows = (block.table.rows || []).some(row =>
            row["FDR"] !== null && row["FDR"] !== "" && isFinite(Number(row["FDR"])) && Number(row["FDR"]) < alpha
        );
        if (!hasSigRows) { nav.finish(); return; }

        function makeChartPane(title) {
            const p = nav.pane(title, { group: "Charts" });
            const controls = el("div", { class: "chart-toggle-row" });
            const area = el("div", { class: "chart-pair" });
            p.appendChild(controls);
            p.appendChild(area);
            const toggles = [];
            return {
                el: p,
                option: (label, options, initialValue, onChange) => {
                    const control = el("label", { class: "chart-option-control" });
                    control.appendChild(el("span", { class: "chart-option-label" }, label));
                    const sel = el("select", { class: "chart-option-select" });
                    options.forEach(([optionLabel, optionValue]) => {
                        sel.appendChild(el("option", { value: optionValue }, optionLabel));
                    });
                    sel.value = initialValue;
                    sel.addEventListener("change", () => onChange(sel.value));
                    control.appendChild(sel);
                    controls.appendChild(control);
                    return sel;
                },
                toggle: (kind, label, buildFn) => {
                    const chartId = _oraChartId(kind, slug, null);
                    const isNetwork = kind === "enrichmap" || kind === "network";
                    const chartDiv = el("div", { id: chartId, class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
                    const chartBtn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
                    let chartShown = false;
                    let chartBuilt = false;
                    function renderChart() { buildFn(chartDiv); chartBuilt = true; }
                    chartBtn.addEventListener("click", () => {
                        if (!chartShown) {
                            chartDiv.style.display = "";
                            if (!chartBuilt) renderChart();
                            else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                            chartShown = true;
                            chartBtn.textContent = "Hide " + label;
                        } else {
                            chartDiv.style.display = "none";
                            chartShown = false;
                            chartBtn.textContent = "Show " + label;
                        }
                    });
                    controls.appendChild(chartBtn);
                    area.appendChild(chartDiv);
                    const t = { redraw: () => { if (chartShown) renderChart(); } };
                    toggles.push(t);
                    return t;
                },
                redraw: () => toggles.forEach(t => t.redraw()),
            };
        }

        const overview = makeChartPane("Protein set overview");
        overview.el.insertBefore(el("div", { class: "ora-chart-threshold-note" },
            "Significance threshold: FDR < " + alpha + " (used for table highlighting, chart filtering, and protein-volcano coloring; adjustable above)."),
            overview.el.firstChild);
        let dotColorScale1 = ORA_DOT_COLOR_SCALES[0][1];
        let barColor = ORA_BAR_COLORS[0][1];
        let topN = 20;
        let sortBy = "fdr";
        let fontSize = 12;

        overview.option("Terms shown", [
            ["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"],
            ["50", "50"], ["100", "100"], ["200", "200"], ["All significant", "all"],
        ], String(topN === Infinity ? "all" : topN), (value) => {
            topN = value === "all" ? Infinity : Number(value);
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Rank & show by", [["FDR", "fdr"], ["Fold enrichment", "fold"], ["Protein count", "k"], ["Protein ratio", "ratio"]], sortBy, (value) => {
            sortBy = value;
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], String(fontSize), (value) => {
            fontSize = Number(value);
            overview.redraw();
            relationships.redraw();
        });

        const dotToggle = overview.toggle("dot", "dot plot", (div) => buildDotPlot(div, block.table, slug, null, dotColorScale1, { topN, sortBy, fontSize, alpha }));
        overview.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale1, (value) => {
            dotColorScale1 = value;
            dotToggle.redraw();
        });
        const barToggle = overview.toggle("bar", "bar plot", (div) => buildBarPlot(div, block.table, slug, null, barColor, { topN, sortBy, fontSize, alpha }));
        overview.option("Bar color", ORA_BAR_COLORS, barColor, (value) => {
            barColor = value;
            barToggle.redraw();
        });
        overview.toggle("proteinvolcano", "protein volcano", (div) =>
            buildTranscriptomicsGeneVolcano(div, resolvedRows, slug, padjCutoff, lfcCutoff, { fontSize }));

        const relationships = makeChartPane("Protein set relationships");
        relationships.el.insertBefore(el("div", { class: "ora-chart-threshold-note ora-chart-encoding-note" },
            "Circle color = -log10(FDR) (darker = more significant); circle size = number of significant proteins in the term. " +
            "Rank & show by and Terms shown (Protein set overview) choose which terms are drawn."),
            relationships.el.firstChild);
        let dotColorScale2 = ORA_DOT_COLOR_SCALES[0][1];
        let edgeThreshold = 1;
        const enrichMapToggle = relationships.toggle("enrichmap", "enrichment map", (div) =>
            buildEnrichmentMap(div, block.table, slug, null, dotColorScale2, { topN, sortBy, fontSize, edgeThreshold, alpha }));
        relationships.option("Edge threshold (shared proteins)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (value) => {
            edgeThreshold = Number(value);
            enrichMapToggle.redraw();
        });
        const networkToggle = relationships.toggle("network", "protein-concept network", (div) =>
            buildTranscriptomicsGeneConceptNetwork(div, block.table, slug, proteinLog2FC, { topN, sortBy, fontSize, alpha, dotColorScale: dotColorScale2 }));
        relationships.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale2, (value) => {
            dotColorScale2 = value;
            enrichMapToggle.redraw();
            networkToggle.redraw();
        });

        nav.finish();
    }

    select.addEventListener("change", () => renderPanel(select.value));
    if (sortedKeys.length > 0) {
        select.value = sortedKeys[0];
        renderPanel(sortedKeys[0]);
    }

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.proteomics.oraConfig = null;
        STATE.proteomics.oraResults = null;
        STATE.proteomics.oraRunLocked = false;
        STATE.proteomics._oraDraft = null;
        STATE.aaPage = "proteomics_ora_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToProteomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// Transcriptomics — Gene Set Enrichment Analysis (GSEA, pre-ranked)
//
// Retrieval (upload/resolve) and the custom-database/background-list widgets
// are FULLY SHARED with ORA above (renderTranscriptomicsCustomDbWidget /
// renderTranscriptomicsBackgroundWidget hit the same server-side session
// keys and are reused here unchanged). Bidirectional by construction - NES
// sign IS the direction - so unlike ORA there is no all/up/down split and
// result keys carry no direction suffix ("GO: Biological Process", not
// "GO: Biological Process (all)").
// ---------------------------------------------------------------------------

const GSEA_RANKING_METRIC_LABELS = {
    signed_log_pvalue: "Signed -log10(p-value)",
    log2fc: "log2 Fold-Change",
    pi_value: "π-value",
    signed_log_padj: "Signed -log10(padj)",
};
// GSEA's own convention (FDR q < 0.25) is looser than ORA's 0.05, since GSEA's
// permutation-pooled FDR is already a stricter correction on its own.
const GSEA_DEFAULT_SIG_THRESHOLD = 0.25;

const TRANSCRIPTOMICS_GSEA_FAMILY_ORDER = [
    "GO: Biological Process", "GO: Molecular Function", "GO: Cellular Component",
    "Reactome",
    ...TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.map((c) => "MSigDB: " + c.label),
    "Custom database",
];
function _sortedTranscriptomicsGseaKeys(tables) {
    return Object.keys(tables).sort((a, b) => {
        const iA = TRANSCRIPTOMICS_GSEA_FAMILY_ORDER.indexOf(a);
        const iB = TRANSCRIPTOMICS_GSEA_FAMILY_ORDER.indexOf(b);
        return (iA === -1 ? 999 : iA) - (iB === -1 ? 999 : iB);
    });
}
function _transcriptomicsGseaSlug(key) {
    return key.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "");
}

function renderTranscriptomicsGseaConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Gene Set Enrichment Analysis: setup"));

    const locked = STATE.transcriptomics.gseaRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    if (!STATE.transcriptomics._gseaDraft) {
        STATE.transcriptomics._gseaDraft = {
            databases: [],
            msigdbCollections: [],
            metric: "signed_log_pvalue",
            minTermSize: 15,
            maxTermSize: 500,
            nPermutations: 1000,
            weight: 1.0,
            seed: 42,
            includeIea: true,
            useBackground: false,
        };
    }
    const draft = STATE.transcriptomics._gseaDraft;

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const runBtnWrap = el("div", { class: "col" });

    function renderGseaRunBtn() {
        runBtnWrap.innerHTML = "";
        const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
        const needsCustomDb = draft.databases.indexOf("Custom database") !== -1 && !STATE.transcriptomics.customDb;
        const needsMsigdbCollection = draft.databases.indexOf("MSigDB") !== -1 && draft.msigdbCollections.length === 0;
        runBtn.disabled = locked || draft.databases.length === 0 || needsCustomDb || needsMsigdbCollection;
        runBtn.addEventListener("click", () => startTranscriptomicsGseaRun(draft));
        runBtnWrap.appendChild(runBtn);
    }

    // Same 4 databases and same MSigDB human/mouse-only gate as ORA
    // (TRANSCRIPTOMICS_ORA_DATABASES, defined above in the ORA section).
    const taxonId = STATE.transcriptomics.taxonId;
    const msigdbAvailable = _msigdbAvailableForTaxon(taxonId);
    TRANSCRIPTOMICS_ORA_DATABASES.forEach((dbName) => {
        const unavailable = dbName === "MSigDB" && !msigdbAvailable;
        if (unavailable) {
            draft.databases = draft.databases.filter((d) => d !== dbName);
        }
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked || unavailable;
        cb.checked = !unavailable && draft.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.databases.indexOf(dbName) === -1) draft.databases.push(dbName);
            } else {
                draft.databases = draft.databases.filter((d) => d !== dbName);
            }
            if (dbName === "Custom database" || dbName === "MSigDB") render();
            else renderGseaRunBtn();
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName + (unavailable ? " (human/mouse only)" : "")));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    if (!msigdbAvailable) {
        wrap.appendChild(el("div", { class: "caption" },
            "MSigDB is disabled: it covers human and mouse only. GO and Reactome support this organism."));
    }
    if (draft.databases.indexOf("MSigDB") !== -1 && msigdbAvailable) {
        wrap.appendChild(renderMsigdbCollectionSelector(locked, draft, taxonId, renderGseaRunBtn));
    }
    if (draft.databases.indexOf("Custom database") !== -1) {
        wrap.appendChild(renderTranscriptomicsCustomDbWidget(locked));
    }
    wrap.appendChild(renderTranscriptomicsBackgroundWidget(locked, draft));

    wrap.appendChild(el("div", { class: "section-label" }, "Ranking metric"));
    const metricRow = el("div", { class: "row" });
    const metricCol = el("div", { class: "col" });
    const metricSelect = el("select", { class: "text-input" });
    [
        ["Signed -log10(p-value) (recommended)", "signed_log_pvalue"],
        ["log2 Fold-Change", "log2fc"],
        ["π-value (log2FC × -log10 p)", "pi_value"],
        ["Signed -log10(padj)", "signed_log_padj"],
    ].forEach(([lbl, val]) => metricSelect.appendChild(el("option", { value: val }, lbl)));
    metricSelect.value = draft.metric;
    metricSelect.disabled = locked;
    metricSelect.addEventListener("change", () => { draft.metric = metricSelect.value; });
    metricCol.appendChild(metricSelect);
    metricRow.appendChild(metricCol);
    wrap.appendChild(metricRow);
    wrap.appendChild(el("div", { class: "caption" },
        "Every gene is ranked by this metric, then GSEA looks for gene sets clustered near either end of " +
        "the ranking. Signed -log10(p-value) (sign of log2FC × significance) is the standard choice " +
        "for a pre-ranked run from a differential-expression table."));

    function numberField(labelText, value, step, minVal, onChange) {
        const col = el("div", { class: "col" });
        col.appendChild(el("div", { class: "chart-option-label" }, labelText));
        const input = el("input", { type: "number", class: "text-input", step: step, min: minVal });
        input.value = value;
        input.disabled = locked;
        input.addEventListener("change", () => {
            const v = parseFloat(input.value);
            if (isFinite(v)) onChange(v);
        });
        col.appendChild(input);
        return col;
    }

    const termSizeRow = el("div", { class: "row" });
    const minCol = el("div", { class: "col" });
    minCol.appendChild(el("div", { class: "chart-option-label" }, "Minimum gene-set size"));
    const minInput = el("input", { type: "number", class: "text-input", min: "2", step: "1" });
    minInput.value = draft.minTermSize;
    minInput.disabled = locked;
    minInput.addEventListener("change", () => {
        const v = parseInt(minInput.value, 10);
        draft.minTermSize = (isFinite(v) && v >= 2) ? v : 15;
        minInput.value = draft.minTermSize;
    });
    minCol.appendChild(minInput);
    termSizeRow.appendChild(minCol);

    const maxCol = el("div", { class: "col" });
    maxCol.appendChild(el("div", { class: "chart-option-label" }, "Maximum gene-set size (blank = no limit)"));
    const maxInput = el("input", { type: "number", class: "text-input", min: "1", step: "1", placeholder: "no limit" });
    maxInput.value = draft.maxTermSize;
    maxInput.disabled = locked;
    maxInput.addEventListener("change", () => {
        const raw = maxInput.value.trim();
        draft.maxTermSize = raw === "" ? "" : parseInt(raw, 10);
    });
    maxCol.appendChild(maxInput);
    termSizeRow.appendChild(maxCol);
    wrap.appendChild(termSizeRow);
    wrap.appendChild(el("div", { class: "caption" },
        "GSEA's own convention (15/500 - stricter at the low end than ORA's 10) since enrichment scores are " +
        "unstable for very small gene sets."));

    const paramRow = el("div", { class: "row" });
    const permCol = el("div", { class: "col" });
    permCol.appendChild(el("div", { class: "chart-option-label" }, "Number of permutations"));
    const permSelect = el("select", { class: "text-input" });
    [["500", "500"], ["1000", "1000"], ["2000", "2000"], ["5000", "5000"]].forEach(([lbl, val]) =>
        permSelect.appendChild(el("option", { value: val }, lbl)));
    permSelect.value = String(draft.nPermutations);
    permSelect.disabled = locked;
    permSelect.addEventListener("change", () => { draft.nPermutations = Number(permSelect.value); });
    permCol.appendChild(permSelect);
    paramRow.appendChild(permCol);
    paramRow.appendChild(numberField("Weight exponent (p)", draft.weight, "0.1", "0", (v) => { draft.weight = v; }));
    paramRow.appendChild(numberField("Random seed (reproducibility)", draft.seed, "1", null, (v) => { draft.seed = Math.round(v); }));
    wrap.appendChild(paramRow);
    wrap.appendChild(el("div", { class: "caption" },
        "More permutations give a finer p-value floor (1 / (permutations + 1)) at the cost of a longer run. " +
        "Weight 1.0 is the standard GSEA weighting (genes contribute proportional to their ranking-metric " +
        "magnitude). The seed makes a re-run with identical settings reproduce identical results exactly."));

    const ieaRow = el("div", { class: "row", style: { alignItems: "center", marginTop: "0.5rem" } });
    const ieaLabel = el("label", { style: { display: "flex", alignItems: "center", gap: "0.4rem" } });
    const ieaCheckbox = el("input", { type: "checkbox", class: "transcriptomics-gsea-iea-checkbox" });
    ieaCheckbox.checked = draft.includeIea;
    ieaCheckbox.disabled = locked;
    ieaCheckbox.addEventListener("change", () => { draft.includeIea = ieaCheckbox.checked; });
    ieaLabel.appendChild(ieaCheckbox);
    ieaLabel.appendChild(document.createTextNode(" Include electronic annotations (IEA)"));
    ieaRow.appendChild(ieaLabel);
    wrap.appendChild(ieaRow);
    wrap.appendChild(el("div", { class: "caption" },
        "Uncheck to restrict GO and Reactome to curated evidence only. No effect on MSigDB. " +
        "Reactome's curated (non-IEA) coverage is essentially human-only."));

    const actions = el("div", { class: "right-actions", style: { marginTop: "1rem" } });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "transcriptomics_module"; render(); });
    btnRow.appendChild(backBtn);
    renderGseaRunBtn();
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function startTranscriptomicsGseaRun(draft) {
    if (draft.databases.length === 0) return;
    const payload = {
        databases: draft.databases,
        // See the ORA run's comment for why this is gated on selection.
        msigdb_collections: draft.databases.indexOf("MSigDB") !== -1 ? draft.msigdbCollections : [],
        metric: draft.metric,
        min_term_size: draft.minTermSize,
        max_term_size: draft.maxTermSize === "" ? null : draft.maxTermSize,
        n_permutations: draft.nPermutations,
        weight: draft.weight,
        seed: draft.seed,
        include_iea: draft.includeIea,
        use_background: !!(draft.useBackground && STATE.transcriptomics.background),
    };
    apiPostJson("/api/transcriptomics/gsea/configure", payload).then(() => {
        STATE.transcriptomics.gseaConfig = payload;
        STATE.transcriptomics.gseaRunLocked = true;
        STATE.transcriptomics.gseaResults = null;
        STATE.aaPage = "transcriptomics_gsea_results";
        render();
        return apiPostJson("/api/transcriptomics/gsea/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("transcriptomics-gsea-progress-inner");
                const label = document.getElementById("transcriptomics-gsea-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => {
                STATE.transcriptomics.gseaResults = result;
                render();
            },
            onError: (err) => {
                alert("GSEA run failed: " + err);
                STATE.transcriptomics.gseaRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.transcriptomics.gseaRunLocked = false;
        render();
    });
}

function renderTranscriptomicsGseaResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Gene Set Enrichment Analysis: results"));

    if (!STATE.transcriptomics.gseaResults) {
        wrap.appendChild(alertBox("info", "Running analysis... (permutation testing can take from a few seconds to a couple of minutes, depending on how many databases are selected)"));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "transcriptomics-gsea-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "transcriptomics-gsea-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const { tables, summary, ranking_profile, ranking_full } = STATE.transcriptomics.gseaResults;

    wrap.appendChild(metricCards([
        ["Genes ranked", summary.n_genes_ranked],
        ["Ranking metric", GSEA_RANKING_METRIC_LABELS[summary.metric] || summary.metric],
        ["Permutations", summary.n_permutations],
        ["Ties in ranking", summary.n_ties],
    ], 4));
    if (summary.n_dropped_padj > 0) {
        wrap.appendChild(alertBox("info",
            summary.n_dropped_padj + " gene(s) had no adjusted p-value (padj) and were excluded from the " +
            "ranking (only relevant when the \"Signed -log10(padj)\" metric is selected)."));
    }
    // See the identical hook on the ORA results screen for why this exists.
    let refreshMsigdbSummary = null;

    const alphaRow = el("div", { class: "ora-sig-threshold-row" });
    alphaRow.appendChild(el("span", { class: "chart-option-label" }, "Significance threshold (FDR q-value <)"));
    const alphaInput = el("input", { type: "number", class: "ora-sig-threshold-input", min: "0", max: "1", step: "0.01" });
    alphaInput.value = STATE.transcriptomics.gseaSigAlpha;
    alphaInput.addEventListener("change", () => {
        let v = parseFloat(alphaInput.value);
        if (!isFinite(v) || v <= 0 || v > 1) v = GSEA_DEFAULT_SIG_THRESHOLD;
        STATE.transcriptomics.gseaSigAlpha = v;
        alphaInput.value = v;
        renderPanel(select.value);
        if (refreshMsigdbSummary) refreshMsigdbSummary();
    });
    alphaRow.appendChild(alphaInput);
    wrap.appendChild(alphaRow);
    wrap.appendChild(el("div", { class: "caption" },
        "GSEA convention default is FDR q < 0.25 (looser than ORA's 0.05), since GSEA's own permutation " +
        "FDR is already a more conservative correction; adjust freely."));

    // Same multiple-testing caveat as ORA's results screen - see the comment
    // there for the full rationale.
    if (typeof summary.n_independent_families === "number") {
        wrap.appendChild(el("div", { class: "caption" },
            "This run scanned " + summary.n_independent_families + " independently FDR-corrected " +
            "result table(s) (no correction is applied across them). Selecting more databases/MSigDB " +
            "collections and reporting the best hit across all of them inflates the true error rate " +
            "above what any single FDR column states."));
    }

    // Cross-collection summary: which MSigDB collection is carrying the
    // signal, at a glance, without clicking through every rail entry.
    const hasMsigdbFamilies = Object.keys(tables).some((k) => k.startsWith("MSigDB: "));
    if (hasMsigdbFamilies) {
        // Card-wrapped - see the identical comment on the ORA results screen.
        const summaryCard = el("div", { class: "ora-db-card msigdb-summary-card" });
        summaryCard.appendChild(el("div", { class: "result-pane-header" },
            el("h2", { class: "result-pane-title" }, "MSigDB collections")));
        const summaryDiv = el("div", { id: _oraChartId("msigdbsummary", "transcriptomics_gsea", null), class: "chart-container", style: { display: "none" } });
        let summaryShown = false, summaryBuilt = false;
        const summaryBtn = el("button", { class: "btn chart-toggle-btn" }, "Show MSigDB collection summary");
        summaryBtn.addEventListener("click", () => {
            if (!summaryShown) {
                summaryDiv.style.display = "";
                if (!summaryBuilt) {
                    buildMsigdbGseaCollectionSummary(summaryDiv, tables, STATE.transcriptomics.gseaSigAlpha, "transcriptomics_gsea");
                    summaryBuilt = true;
                } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(summaryDiv);
                summaryShown = true; summaryBtn.textContent = "Hide MSigDB collection summary";
            } else { summaryDiv.style.display = "none"; summaryShown = false; summaryBtn.textContent = "Show MSigDB collection summary"; }
        });
        refreshMsigdbSummary = () => {
            summaryBuilt = false;
            if (summaryShown) {
                buildMsigdbGseaCollectionSummary(summaryDiv, tables, STATE.transcriptomics.gseaSigAlpha, "transcriptomics_gsea");
                summaryBuilt = true;
            }
        };
        summaryCard.appendChild(summaryBtn);
        summaryCard.appendChild(summaryDiv);
        wrap.appendChild(summaryCard);
    }

    const sortedKeys = _sortedTranscriptomicsGseaKeys(tables);
    const selectorRow = el("div", { class: "row" });
    const select = el("select", { class: "text-input transcriptomics-gsea-result-select" });
    sortedKeys.forEach((key) => {
        const block = tables[key];
        const suffix = block.n_terms_tested > 0 ? ` (${block.n_terms_tested} gene sets)` : " (no gene sets)";
        select.appendChild(el("option", { value: key }, key + suffix));
    });
    selectorRow.appendChild(el("div", { class: "col" }, select));
    wrap.appendChild(selectorRow);

    const panelHost = el("div", {});
    wrap.appendChild(panelHost);

    function renderPanel(key) {
        panelHost.innerHTML = "";
        const block = tables[key];
        if (!block) return;
        const alpha = STATE.transcriptomics.gseaSigAlpha;
        const slug = _transcriptomicsGseaSlug(key);
        const card = el("div", { class: "rank-card ora-db-card ora-result-card ora-db-card-" + slug });
        panelHost.appendChild(card);
        const navHost = el("div", {});
        card.appendChild(navHost);
        const nav = makeResultNav(navHost, { stateKey: "transcriptomics_gsea" });

        const overviewPane = nav.pane("Overview");
        if (block.note) overviewPane.appendChild(alertBox("info", block.note));
        if (block.n_terms_tested > 0) {
            overviewPane.appendChild(metricCards([
                ["Gene sets tested", block.n_terms_tested],
                ["Universe (ranked genes)", block.universe_size],
                ["Excluded (size filter)", block.n_terms_excluded_min_size],
            ], 3));
        }
        if (!block.table.rows || block.table.rows.length === 0) {
            if (!block.note) {
                overviewPane.appendChild(alertBox("info", "No gene sets passed the size filter for this database."));
            }
            nav.finish();
            return;
        }

        const tablePane = nav.pane("Result table");
        tablePane.appendChild(buildTable(block.table, {
            maxHeight: 480,
            rowClassFn: (row) => {
                const fdr = row["FDR q-value"];
                return (fdr !== "" && fdr !== null && isFinite(Number(fdr)) && Number(fdr) < alpha) ? "row-significant" : null;
            },
            cellFormatFn: (col, val) => {
                if (col === "p-value" || col === "FDR q-value" || col === "FDR (BH)") return _fmtSci(val, 3);
                if (col === "ES" || col === "NES" || col === "Leading Edge Fraction") return _fmtFixed(val, 3);
                return val;
            },
        }));
        const dlBtn = el("button", { class: "btn" }, "Download result table");
        dlBtn.addEventListener("click", () => {
            // Per-key filename - see the ORA download's comment for why a
            // fixed name silently collides across every result family.
            downloadFile("/api/transcriptomics/gsea/download?key=" + encodeURIComponent(key),
                         "transcriptomics_gsea_" + _transcriptomicsGseaSlug(key) + ".csv");
        });
        tablePane.appendChild(dlBtn);

        const hasSigRows = (block.table.rows || []).some(row =>
            row["FDR q-value"] !== null && row["FDR q-value"] !== "" && isFinite(Number(row["FDR q-value"])) && Number(row["FDR q-value"]) < alpha
        );
        if (!hasSigRows) { nav.finish(); return; }

        // Generic container-scoped option/toggle factories (each group builds
        // its own row/area pair, so a control only ever redraws its own
        // group's charts - a shared flat row was the problem this replaced).
        function addChartOption(rowEl, label, options, initialValue, onChange) {
            const control = el("label", { class: "chart-option-control" });
            control.appendChild(el("span", { class: "chart-option-label" }, label));
            const sel = el("select", { class: "chart-option-select" });
            options.forEach(([optionLabel, optionValue]) => {
                sel.appendChild(el("option", { value: optionValue }, optionLabel));
            });
            sel.value = initialValue;
            sel.addEventListener("change", () => onChange(sel.value));
            control.appendChild(sel);
            rowEl.appendChild(control);
            return sel;
        }

        function addChartToggle(rowEl, areaEl, kind, label, buildFn) {
            const chartId = _oraChartId(kind, slug, null);
            const isNetwork = kind === "enrichmap" || kind === "network";
            const chartDiv = el("div", { id: chartId, class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
            const chartBtn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
            let chartShown = false;
            let chartBuilt = false;
            function renderChart() { buildFn(chartDiv); chartBuilt = true; }
            chartBtn.addEventListener("click", () => {
                if (!chartShown) {
                    chartDiv.style.display = "";
                    if (!chartBuilt) renderChart();
                    else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                    chartShown = true;
                    chartBtn.textContent = "Hide " + label;
                } else {
                    chartDiv.style.display = "none";
                    chartShown = false;
                    chartBtn.textContent = "Show " + label;
                }
            });
            rowEl.appendChild(chartBtn);
            areaEl.appendChild(chartDiv);
            return { redraw: () => { if (chartShown) renderChart(); } };
        }

        // Each chart group is now its own rail entry nested under "Charts",
        // rather than a stacked section on one long scroll.
        function makeChartGroup(title, caption) {
            const p = nav.pane(title, { group: "Charts" });
            if (caption) p.appendChild(el("div", { class: "caption" }, caption));
            const controls = el("div", { class: "chart-toggle-row" });
            const area = el("div", { class: "chart-pair" });
            p.appendChild(controls);
            p.appendChild(area);
            const toggles = [];
            return {
                el: p,
                option: (l, o, i, f) => addChartOption(controls, l, o, i, f),
                toggle: (k, l, b) => { const t = addChartToggle(controls, area, k, l, b); toggles.push(t); return t; },
                redraw: () => toggles.forEach(t => t.redraw()),
            };
        }

        // Shared, cross-group state (genuinely affects every figure) - sits
        // above the section rail so it stays visible no matter which section
        // is open, rather than living inside any one pane.
        let topN = 20;
        let fontSize = 12;
        const groupsAffectedByTopN = [];   // overview + relationships + leadingEdge (NOT per-gene-set or diagnostics)
        const allGroupsForFontSize = [];   // every group

        const displayRow = el("div", { class: "chart-toggle-row result-display-row" });
        card.insertBefore(displayRow, navHost);
        addChartOption(displayRow, "Gene sets shown", [
            ["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"],
            ["50", "50"], ["All significant", "all"],
        ], String(topN === Infinity ? "all" : topN), (value) => {
            topN = value === "all" ? Infinity : Number(value);
            groupsAffectedByTopN.forEach(g => g.redraw());
        });
        addChartOption(displayRow, "Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], String(fontSize), (value) => {
            fontSize = Number(value);
            allGroupsForFontSize.forEach(g => g.redraw());
        });

        // --- Group 1: Gene set overview -------------------------------------
        const overview = makeChartGroup("Gene set overview");
        overview.el.insertBefore(el("div", { class: "ora-chart-threshold-note" },
            "Significance threshold: FDR q-value < " + alpha + " (used for table highlighting, chart filtering, and volcano coloring)."),
            overview.el.firstChild);
        let dotColorScale1 = ORA_DOT_COLOR_SCALES[0][1];
        const dotToggle = overview.toggle("dot", "NES dot plot", (div) => buildGseaDotPlot(div, block.table, slug, dotColorScale1, { topN, fontSize, alpha }));
        overview.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale1, (value) => {
            dotColorScale1 = value;
            dotToggle.redraw();
        });
        overview.toggle("bar", "NES bar plot", (div) => buildGseaBarPlot(div, block.table, slug, { topN, fontSize, alpha }));
        overview.toggle("volcano", "NES volcano plot", (div) => buildGseaVolcanoPlot(div, block.table, slug, alpha, { fontSize }));
        groupsAffectedByTopN.push(overview);
        allGroupsForFontSize.push(overview);

        // --- Group 2: Gene set relationships ---------------------------------
        const relationships = makeChartGroup("Gene set relationships");
        let dotColorScale2 = ORA_DOT_COLOR_SCALES[0][1];
        let edgeThreshold = 1;
        const enrichMapToggle = relationships.toggle("enrichmap", "enrichment map", (div) =>
            buildGseaEnrichmentMap(div, block.table, slug, dotColorScale2, { topN, fontSize, edgeThreshold, alpha }));
        relationships.option("Edge threshold (shared genes)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (value) => {
            edgeThreshold = Number(value);
            enrichMapToggle.redraw();
        });
        const networkToggle = relationships.toggle("network", "leading-edge gene network", (div) =>
            buildGseaGeneConceptNetwork(div, block.table, slug, { topN, fontSize, alpha, dotColorScale: dotColorScale2 }));
        relationships.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale2, (value) => {
            dotColorScale2 = value;
            enrichMapToggle.redraw();
            networkToggle.redraw();
        });
        groupsAffectedByTopN.push(relationships);
        allGroupsForFontSize.push(relationships);

        // --- Group 3: Per gene set (enrichment plot, null histogram) --------
        // Both reuse ONE stable container div across selection changes -
        // mirrors netprop's per-module enrichment pattern (_netpropChartToggle
        // with an explicitId), so both go through it rather than the group
        // factory above (which always creates a FRESH div per toggle).
        const curveTermIds = Object.keys(block.curves || {});
        if (curveTermIds.length > 0) {
            const rowsByTermId = {};
            (block.table.rows || []).forEach(r => { rowsByTermId[r["Term ID"]] = r; });

            const perSetPane = nav.pane("Per gene set", { group: "Charts" });
            perSetPane.appendChild(el("div", { class: "caption" },
                "These two figures show ONE gene set at a time - the one picked in the \"Gene set\" menu " +
                "directly below. Only the top " + curveTermIds.length + " gene set(s) by FDR q-value are " +
                "available here."));

            let selectedTermId = curveTermIds[0];
            let perSetFontSize = 12;
            const perSetControls = el("div", { class: "chart-toggle-row" });
            perSetPane.appendChild(perSetControls);

            const termSelect = el("select", { class: "chart-option-select" });
            curveTermIds.forEach((tid) => {
                const r = rowsByTermId[tid];
                const nesLabel = r ? " (NES " + Number(r.NES).toFixed(2) + ")" : "";
                // Term NAME, falling back to the ID - a GO accession like
                // GO:0006281 is unreadable in a picker, and the figure titles
                // these drive are labelled by name, so the two must agree.
                const name = r && r["Term Name"] ? String(r["Term Name"]) : tid;
                termSelect.appendChild(el("option", { value: tid }, name + nesLabel));
            });
            termSelect.value = selectedTermId;
            termSelect.addEventListener("change", () => {
                selectedTermId = termSelect.value;
                enrichToggle.redraw();
                nullHistToggle.redraw();
            });
            const termControl = el("label", { class: "chart-option-control" });
            termControl.appendChild(el("span", { class: "chart-option-label" }, "Gene set (drives both figures below)"));
            termControl.appendChild(termSelect);
            perSetControls.appendChild(termControl);

            addChartOption(perSetControls, "Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], "12", (value) => {
                perSetFontSize = Number(value);
                enrichToggle.redraw();
                nullHistToggle.redraw();
            });

            const perSetArea = el("div", { class: "chart-pair" });
            perSetPane.appendChild(perSetArea);

            const enrichToggle = _netpropChartToggle("Show enrichment plot", "enrichplot",
                (div) => buildGseaEnrichmentPlot(div, block.curves[selectedTermId], ranking_profile, rowsByTermId[selectedTermId], perSetFontSize, slug),
                _oraChartId("enrichplot", slug, null));
            perSetArea.appendChild(enrichToggle.el);

            const nullHistToggle = _netpropChartToggle("Show null distribution", "nullhist",
                (div) => buildGseaNullHistogram(div, block, rowsByTermId[selectedTermId], perSetFontSize, slug),
                _oraChartId("nullhist", slug, null));
            perSetArea.appendChild(nullHistToggle.el);

            // --- Compare gene sets (the overlay, and ONLY the overlay) - its
            // own rail entry, so its checkboxes can no longer read as if they
            // control the two figures above.
            // Multi-select via checkbox pills (already in the design system via
            // .ora-db-pill), not <select multiple> - a long GO term list is
            // unusable with ctrl-click discovery. Capped at 5 curves: beyond
            // that the running-ES lines are no longer separable by eye.
            const comparePane = nav.pane("Compare gene sets", { group: "Charts" });
            comparePane.appendChild(el("div", { class: "caption" },
                "Overlays several gene sets' enrichment curves on one figure. Tick at least 2."));

            const overlaySelectedIds = [];
            const overlayCountLabel = el("div", { class: "caption" }, "");
            function updateOverlayCount() {
                const n = overlaySelectedIds.length;
                overlayCountLabel.textContent = n >= 5
                    ? "5 of 5 selected - the maximum. Untick one to swap in a different gene set."
                    : n + " of 5 selected" + (n < 2 ? " (tick at least 2 to draw the figure)." : ".");
            }
            const overlayPillsWrap = el("div", { class: "db-checkbox-row" });
            curveTermIds.forEach((tid) => {
                const r = rowsByTermId[tid];
                const label = el("label", { class: "ora-db-pill" });
                const cb = el("input", { type: "checkbox" });
                cb.addEventListener("change", () => {
                    if (cb.checked) {
                        if (overlaySelectedIds.length >= 5) { cb.checked = false; return; }
                        overlaySelectedIds.push(tid);
                    } else {
                        const i = overlaySelectedIds.indexOf(tid);
                        if (i !== -1) overlaySelectedIds.splice(i, 1);
                    }
                    overlayPillsWrap.querySelectorAll("input[type=checkbox]").forEach((el2) => {
                        if (!el2.checked) el2.disabled = overlaySelectedIds.length >= 5;
                    });
                    updateOverlayCount();
                    multiEnrichToggle.redraw();
                });
                label.appendChild(cb);
                label.appendChild(document.createTextNode(
                    " " + (r && r["Term Name"] ? String(r["Term Name"]) : tid) +
                    (r ? " (NES " + Number(r.NES).toFixed(2) + ")" : "")));
                overlayPillsWrap.appendChild(label);
            });
            updateOverlayCount();
            comparePane.appendChild(overlayCountLabel);
            comparePane.appendChild(overlayPillsWrap);

            const overlayArea = el("div", { class: "chart-pair" });
            comparePane.appendChild(overlayArea);
            const multiEnrichToggle = _netpropChartToggle("Show multi-set enrichment plot", "multienrich",
                (div) => buildGseaMultiEnrichmentPlot(div, block.curves, ranking_profile, rowsByTermId, overlaySelectedIds, perSetFontSize, slug),
                _oraChartId("multienrich", slug, null));
            overlayArea.appendChild(multiEnrichToggle.el);
        }

        // --- Group 4: Leading edge (ridgeline, heatmap) ----------------------
        if (ranking_full && ranking_full.n_genes > 0) {
            const leadingEdge = makeChartGroup("Leading edge",
                "Each control below is tagged with the figure it changes - a control tagged \"heatmap\" " +
                "will not move the ridgeline.");
            let population = "edge";
            let geneMode = "count";
            let rowOrder = "nes";
            let maxGenes = 60;
            const ridgeToggle = leadingEdge.toggle("ridgeline", "ridgeline plot", (div) =>
                buildGseaRidgePlot(div, block.table, block.curves, ranking_full, slug, { topN, fontSize, alpha, population }));
            leadingEdge.option("Gene population (ridgeline)", [["Leading edge (core enriched)", "edge"], ["All set members", "members"]], population, (value) => {
                population = value;
                ridgeToggle.redraw();
            });
            const heatmapToggle = leadingEdge.toggle("heatplot", "leading-edge heatmap", (div) =>
                buildGseaLeadingEdgeHeatmap(div, block.table, ranking_full, slug, { topN, fontSize, alpha, geneMode, rowOrder, maxGenes }));
            leadingEdge.option("Genes shown (heatmap)", [["30", "30"], ["45", "45"], ["60", "60"], ["90", "90"]], String(maxGenes), (value) => {
                maxGenes = Number(value);
                heatmapToggle.redraw();
            });
            leadingEdge.option("Gene selection (heatmap)", [["Shared across sets", "count"], ["Strongest metric", "metric"]], geneMode, (value) => {
                geneMode = value;
                heatmapToggle.redraw();
            });
            leadingEdge.option("Row order (heatmap)", [["By NES", "nes"], ["By FDR q", "fdr"]], rowOrder, (value) => {
                rowOrder = value;
                heatmapToggle.redraw();
            });
            groupsAffectedByTopN.push(leadingEdge);
            allGroupsForFontSize.push(leadingEdge);
        }

        // --- Group 5: Diagnostics (QC, not a finding - kept visually separate) ---
        const diagnostics = makeChartGroup("Diagnostics",
            "Quality control, not a biological result.");
        let dotColorScale5 = ORA_DOT_COLOR_SCALES[0][1];
        const nesSizeToggle = diagnostics.toggle("nessize", "NES vs set size", (div) =>
            buildGseaNesVsSizePlot(div, block.table, slug, dotColorScale5, { fontSize }));
        diagnostics.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale5, (value) => {
            dotColorScale5 = value;
            nesSizeToggle.redraw();
        });
        allGroupsForFontSize.push(diagnostics);

        nav.finish();
    }

    select.addEventListener("change", () => renderPanel(select.value));
    if (sortedKeys.length > 0) {
        select.value = sortedKeys[0];
        renderPanel(sortedKeys[0]);
    }

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.transcriptomics.gseaConfig = null;
        STATE.transcriptomics.gseaResults = null;
        STATE.transcriptomics.gseaRunLocked = false;
        STATE.transcriptomics._gseaDraft = null;
        STATE.aaPage = "transcriptomics_gsea_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// Transcriptomics — Transcription Factor / Regulon
//
// Retrieval (upload/resolve) is FULLY SHARED with ORA/GSEA above; unlike
// them, this module defines no custom-database/background widgets of its
// own - the regulon databases themselves ARE the gene-set source (see
// Transcriptomics_regulon.py's module docstring).
//
// Two analyses, either or both selectable:
//   - TARGET ENRICHMENT reuses the exact ORA table shape (Term ID/Term Name/
//     k/K/n/N/Fold Enrichment/p-value/FDR/Bonferroni/Overlapping Genes), so
//     it reuses buildDotPlot/buildBarPlot/buildVolcanoPlot/buildEnrichmentMap/
//     buildTranscriptomicsGeneConceptNetwork UNCHANGED - the whole point of
//     regulons satisfying the same gene-set contract every other database
//     does.
//   - ACTIVITY INFERENCE (TF/Regulon Score column, from the ULM in
//     Transcriptomics_regulon.py) is genuinely new and gets its own two
//     chart builders below (bar plot, volcano).
// ---------------------------------------------------------------------------

function renderTranscriptomicsRegulonConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Transcription Factor / Regulon: setup"));

    const locked = STATE.transcriptomics.regulonRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    if (!STATE.transcriptomics._regulonDraft) {
        STATE.transcriptomics._regulonDraft = {
            databases: [],
            analyses: ["enrichment", "activity"],
            dorotheaLevels: "ABC",
            metric: "signed_log_pvalue",
            padjCutoff: 0.05,
            lfcCutoff: 1.0,
            minTermSize: 5,
            maxTermSize: 500,
            minRegulonSize: 5,
        };
    }
    const draft = STATE.transcriptomics._regulonDraft;
    const taxonId = STATE.transcriptomics.taxonId;

    wrap.appendChild(el("div", { class: "caption" },
        "A regulon is a transcription factor's set of target genes. TRRUST, DoRothEA and CollecTRI record " +
        "whether each target is activated or repressed; ChEA3's libraries do not, so it participates in " +
        "Target Enrichment only."));

    // Warn BEFORE the run, not via an empty results table afterwards. A TF is
    // only testable when enough of its targets are present in the uploaded
    // experiment, so a small resolved table cannot support activity
    // inference at all - a real risk for a genuinely short gene list.
    const nResolvedRows = ((STATE.transcriptomics.resolved
        && STATE.transcriptomics.resolved.resolved_table
        && STATE.transcriptomics.resolved.resolved_table.rows) || []).length;
    if (nResolvedRows > 0 && nResolvedRows < REGULON_MIN_GENES_FOR_USEFUL_RUN) {
        wrap.appendChild(alertBox("warning",
            "Your resolved table has only " + nResolvedRows + " gene(s). A transcription factor is only " +
            "tested when at least \"Minimum regulon size\" of its targets are present in your data, so a " +
            "list this short will very likely produce an EMPTY activity table (and few or no enriched " +
            "TFs). This analysis is designed for a full differential-expression table - typically " +
            "thousands of genes, not a short candidate list. You can still run it, but expect little or " +
            "nothing to come back."));
    }

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const runBtnWrap = el("div", { class: "col" });

    function renderRunBtn() {
        runBtnWrap.innerHTML = "";
        const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
        runBtn.disabled = locked || draft.databases.length === 0 || draft.analyses.length === 0;
        runBtn.addEventListener("click", () => startTranscriptomicsRegulonRun(draft));
        runBtnWrap.appendChild(runBtn);
    }

    let anyGated = false;
    TRANSCRIPTOMICS_REGULON_DATABASES.forEach((dbName) => {
        const available = _regulonDbAvailableForTaxon(dbName, taxonId);
        if (!available) {
            anyGated = true;
            draft.databases = draft.databases.filter((d) => d !== dbName);
        }
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked || !available;
        cb.checked = available && draft.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.databases.indexOf(dbName) === -1) draft.databases.push(dbName);
            } else {
                draft.databases = draft.databases.filter((d) => d !== dbName);
            }
            renderRunBtn();
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName + (available ? "" : " (" + TRANSCRIPTOMICS_REGULON_DB_COVERAGE[dbName] + ")")));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    if (anyGated) {
        wrap.appendChild(el("div", { class: "caption" },
            "Some databases are disabled: TRRUST/DoRothEA/CollecTRI cover human and mouse; ChEA3's libraries " +
            "use human gene symbols and are human only."));
    }

    wrap.appendChild(el("div", { class: "section-label" }, "Analyses to run"));
    const analysisRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    [["enrichment", "Target enrichment"], ["activity", "Activity inference"]].forEach(([key, lbl]) => {
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked;
        cb.checked = draft.analyses.indexOf(key) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (draft.analyses.indexOf(key) === -1) draft.analyses.push(key);
            } else {
                draft.analyses = draft.analyses.filter((a) => a !== key);
            }
            renderRunBtn();
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + lbl));
        analysisRow.appendChild(label);
    });
    wrap.appendChild(analysisRow);

    wrap.appendChild(el("div", { class: "section-label" }, "DoRothEA confidence levels"));
    const dorotheaRow = el("div", { class: "row" });
    const dorotheaSelect = el("select", { class: "text-input" });
    [["A (highest confidence)", "A"], ["A + B", "AB"], ["A + B + C (recommended)", "ABC"], ["A + B + C + D", "ABCD"], ["A + B + C + D + E (all)", "ABCDE"]]
        .forEach(([lbl, val]) => dorotheaSelect.appendChild(el("option", { value: val }, lbl)));
    dorotheaSelect.value = draft.dorotheaLevels;
    dorotheaSelect.disabled = locked;
    dorotheaSelect.addEventListener("change", () => { draft.dorotheaLevels = dorotheaSelect.value; });
    dorotheaRow.appendChild(el("div", { class: "col" }, dorotheaSelect));
    wrap.appendChild(dorotheaRow);

    wrap.appendChild(el("div", { class: "section-label" }, "Ranking metric (activity inference)"));
    const metricRow = el("div", { class: "row" });
    const metricSelect = el("select", { class: "text-input" });
    [
        ["Signed -log10(p-value) (recommended)", "signed_log_pvalue"],
        ["log2 Fold-Change", "log2fc"],
        ["π-value (log2FC × -log10 p)", "pi_value"],
        ["Signed -log10(padj)", "signed_log_padj"],
    ].forEach(([lbl, val]) => metricSelect.appendChild(el("option", { value: val }, lbl)));
    metricSelect.value = draft.metric;
    metricSelect.disabled = locked;
    metricSelect.addEventListener("change", () => { draft.metric = metricSelect.value; });
    metricRow.appendChild(el("div", { class: "col" }, metricSelect));
    wrap.appendChild(metricRow);
    wrap.appendChild(el("div", { class: "caption" },
        "Every gene's movement is measured by this metric before being regressed against each TF's " +
        "activated/repressed targets. Same metric family as GSEA; only used by activity inference."));

    function numberField(labelText, value, step, minVal, onChange) {
        const col = el("div", { class: "col" });
        col.appendChild(el("div", { class: "chart-option-label" }, labelText));
        const input = el("input", { type: "number", class: "text-input", step: step, min: minVal });
        input.value = value;
        input.disabled = locked;
        input.addEventListener("change", () => {
            const v = parseFloat(input.value);
            if (isFinite(v)) onChange(v);
        });
        col.appendChild(input);
        return col;
    }

    wrap.appendChild(el("div", { class: "section-label" }, "Target enrichment thresholds"));
    const degRow = el("div", { class: "row" });
    degRow.appendChild(numberField("padj cutoff", draft.padjCutoff, "0.01", "0", (v) => { draft.padjCutoff = v; }));
    degRow.appendChild(numberField("log2FC cutoff", draft.lfcCutoff, "0.1", "0", (v) => { draft.lfcCutoff = v; }));
    wrap.appendChild(degRow);
    wrap.appendChild(el("div", { class: "caption" },
        "Defines which of your genes count as differentially expressed for target enrichment (Fisher's exact " +
        "test), exactly as in Over-Representation Analysis. Does not affect activity inference, which uses " +
        "every tested gene's continuous ranking metric instead of a DEG cutoff."));

    const sizeRow = el("div", { class: "row" });
    const minTermCol = el("div", { class: "col" });
    minTermCol.appendChild(el("div", { class: "chart-option-label" }, "Minimum regulon size (both analyses)"));
    const minTermInput = el("input", { type: "number", class: "text-input", min: "1", step: "1" });
    minTermInput.value = draft.minTermSize;
    minTermInput.disabled = locked;
    minTermInput.addEventListener("change", () => {
        const v = parseInt(minTermInput.value, 10);
        draft.minTermSize = (isFinite(v) && v >= 1) ? v : 5;
        draft.minRegulonSize = draft.minTermSize;
        minTermInput.value = draft.minTermSize;
    });
    minTermCol.appendChild(minTermInput);
    sizeRow.appendChild(minTermCol);

    const maxTermCol = el("div", { class: "col" });
    maxTermCol.appendChild(el("div", { class: "chart-option-label" }, "Maximum regulon size, enrichment only (blank = no limit)"));
    const maxTermInput = el("input", { type: "number", class: "text-input", min: "1", step: "1", placeholder: "no limit" });
    maxTermInput.value = draft.maxTermSize;
    maxTermInput.disabled = locked;
    maxTermInput.addEventListener("change", () => {
        const raw = maxTermInput.value.trim();
        draft.maxTermSize = raw === "" ? "" : parseInt(raw, 10);
    });
    maxTermCol.appendChild(maxTermInput);
    sizeRow.appendChild(maxTermCol);
    wrap.appendChild(sizeRow);
    wrap.appendChild(el("div", { class: "caption" },
        "A TF with very few targets present in your experiment gives an unreliable result either way - too " +
        "few overlapping genes for enrichment, too few points to fit for activity inference."));

    const actions = el("div", { class: "right-actions", style: { marginTop: "1rem" } });
    const btnRow = el("div", { class: "btn-row" });
    const backBtn = el("button", { class: "btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "transcriptomics_module"; render(); });
    btnRow.appendChild(backBtn);
    renderRunBtn();
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function startTranscriptomicsRegulonRun(draft) {
    if (draft.databases.length === 0 || draft.analyses.length === 0) return;
    const payload = {
        databases: draft.databases,
        analyses: draft.analyses,
        dorothea_levels: draft.dorotheaLevels,
        metric: draft.metric,
        padj_cutoff: draft.padjCutoff,
        lfc_cutoff: draft.lfcCutoff,
        min_term_size: draft.minTermSize,
        max_term_size: draft.maxTermSize === "" ? null : draft.maxTermSize,
        min_regulon_size: draft.minRegulonSize,
    };
    apiPostJson("/api/transcriptomics/regulon/configure", payload).then(() => {
        STATE.transcriptomics.regulonConfig = payload;
        STATE.transcriptomics.regulonRunLocked = true;
        STATE.transcriptomics.regulonResults = null;
        STATE.aaPage = "transcriptomics_regulon_results";
        render();
        return apiPostJson("/api/transcriptomics/regulon/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("transcriptomics-regulon-progress-inner");
                const label = document.getElementById("transcriptomics-regulon-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => {
                STATE.transcriptomics.regulonResults = result;
                render();
            },
            onError: (err) => {
                alert("TF / Regulon run failed: " + err);
                STATE.transcriptomics.regulonRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.transcriptomics.regulonRunLocked = false;
        render();
    });
}

function _transcriptomicsRegulonSlug(key) {
    return "regulon_" + key.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_+|_+$/g, "");
}

function renderTranscriptomicsRegulonResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Transcriptomics — Transcription Factor / Regulon: results"));

    if (!STATE.transcriptomics.regulonResults) {
        wrap.appendChild(alertBox("info", "Running analysis... (fetching regulon databases and testing every transcription factor can take from a few seconds to a couple of minutes)"));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "transcriptomics-regulon-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "transcriptomics-regulon-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const results = STATE.transcriptomics.regulonResults;
    const enrichment = results.enrichment || null;
    const activity = results.activity || null;

    const alphaRow = el("div", { class: "ora-sig-threshold-row" });
    alphaRow.appendChild(el("span", { class: "chart-option-label" }, "Significance threshold (FDR <)"));
    const alphaInput = el("input", { type: "number", class: "ora-sig-threshold-input", min: "0", max: "1", step: "0.001" });
    alphaInput.value = STATE.transcriptomics.regulonSigAlpha;
    alphaInput.addEventListener("change", () => {
        let v = parseFloat(alphaInput.value);
        if (!isFinite(v) || v <= 0 || v > 1) v = 0.05;
        STATE.transcriptomics.regulonSigAlpha = v;
        alphaInput.value = v;
        render();
    });
    alphaRow.appendChild(alphaInput);
    wrap.appendChild(alphaRow);
    // Each result table below is its OWN Benjamini-Hochberg family - the FDR
    // in one table says nothing about any other table. Two databases already
    // produce 6 enrichment tables (database x all/up/down) plus one activity
    // table per database, with no correction across them - a real, easy-to-
    // miss gap found while auditing this module. Same convention ORA uses;
    // stated here because the READMEs are not reachable from the app yet.
    wrap.appendChild(el("div", { class: "caption" },
        "Each result table above is corrected separately (its own Benjamini-Hochberg family) - the " +
        "threshold applies WITHIN a table, not across them. Selecting 2 databases already gives 6 " +
        "enrichment tables (each database x all/up/down) plus one activity table per database; scanning " +
        "all of them and reporting whichever hit looks best inflates your real false-discovery rate " +
        "beyond what any single FDR column shows. Decide which table answers your question before " +
        "looking, not after."));

    // Bordered card with the standard blue top accent, matching GSEA/ORA
    // (app.js's makeResultNav callers) - this screen used to append navHost
    // straight to the page with no box, visually inconsistent with every
    // other results screen. navHost stays a separate nested div (not passed
    // directly to makeResultNav) so nothing above the rail - none needed
    // here, but keeps the pattern identical to GSEA's card.
    const card = el("div", { class: "rank-card ora-db-card ora-result-card" });
    wrap.appendChild(card);
    const navHost = el("div", {});
    card.appendChild(navHost);
    const nav = makeResultNav(navHost, { stateKey: "transcriptomics_regulon" });
    const alpha = STATE.transcriptomics.regulonSigAlpha;

    // Per-gene log2FC lookup (uppercase symbol -> log2FC), built once from the
    // retriever's resolved table - mirrors the identical lookup the ORA
    // results screen builds (app.js, renderTranscriptomicsOraResultsScreen).
    // Feeds the gene-concept network's gene coloring (Figures C/D) and the
    // activity-vs-own-expression scatter (Figure B): a TF's OWN log2FC in
    // this experiment, looked up by its symbol.
    const geneLog2FC = {};
    // ...and the matching per-gene padj, so the activity-vs-own-expression
    // figure can show whether a TF's OWN mRNA change is statistically real,
    // not just how big it is. Without this the y-axis carries magnitude with
    // no confidence attached, which is what made that figure unreadable.
    const genePadj = {};
    const resolvedRowsForRegulon = (STATE.transcriptomics.resolved && STATE.transcriptomics.resolved.resolved_table && STATE.transcriptomics.resolved.resolved_table.rows) || [];
    resolvedRowsForRegulon.forEach((r) => {
        const sym = String(r.gene_symbol || "").toUpperCase();
        if (sym && r.log2FC !== null && isFinite(Number(r.log2FC))) geneLog2FC[sym] = Number(r.log2FC);
        if (sym && r.padj !== null && isFinite(Number(r.padj))) genePadj[sym] = Number(r.padj);
    });

    // --- Pane: TF Activity ---------------------------------------------------
    if (activity) {
        const activityPane = nav.pane("TF Activity");
        const activityDbs = Object.keys(activity.tables);
        if (activityDbs.length === 0) {
            activityPane.appendChild(alertBox("info", "Activity inference was not run for this analysis."));
        } else {
            const selRow = el("div", { class: "row" });
            const dbSelect = el("select", { class: "text-input" });
            activityDbs.forEach((db) => {
                const block = activity.tables[db];
                const suffix = block.n_tfs_tested > 0 ? ` (${block.n_tfs_tested} TFs)` : " (no TFs)";
                dbSelect.appendChild(el("option", { value: db }, db + suffix));
            });
            selRow.appendChild(el("div", { class: "col" }, dbSelect));
            activityPane.appendChild(selRow);
            const activityHost = el("div", {});
            activityPane.appendChild(activityHost);

            function renderActivityDb(db) {
                activityHost.innerHTML = "";
                const block = activity.tables[db];
                if (!block) return;
                if (block.note) activityHost.appendChild(alertBox("info", block.note));
                if (!block.table.rows || block.table.rows.length === 0) { return; }

                const slug = _transcriptomicsRegulonSlug("activity_" + db);
                const actSummary = activity.summary || {};
                activityHost.appendChild(metricCards([
                    ["TFs tested", block.n_tfs_tested],
                    ["Ranking metric", GSEA_RANKING_METRIC_LABELS[actSummary.metric] || actSummary.metric || "-"],
                    ["Permutations", actSummary.n_permutations || "-"],
                    ["Random seed", actSummary.seed === undefined ? "-" : actSummary.seed],
                    ["Null inflation (λ)", Number(block.lambda_used).toFixed(3)],
                ], 5));
                if (block.calibration_note) {
                    activityHost.appendChild(alertBox("info", block.calibration_note));
                }
                activityHost.appendChild(el("div", { class: "caption" },
                    "λ (\"null inflation\") is how much wider this run's spread of TF scores is than plain " +
                    "theory predicts. \"FDR (calibrated)\" divides that factor back out and is what sorts and " +
                    "highlights this table; it can only ever make results MORE conservative, never less."));
                activityHost.appendChild(alertBox("info",
                    "λ is still ambiguous, though less than a plain median/MAD estimate would be. A large λ " +
                    "can mean (a) the null really is inflated because your TFs' targets are co-regulated - " +
                    "in which case calibrating is right - OR (b) a large fraction of your TFs are genuinely " +
                    "active, which widens the spread too, in which case calibration is more conservative " +
                    "than it needs to be. This app uses a truncated-normal estimate rather than plain MAD " +
                    "specifically because MAD is fooled by (b) badly (measured: 25% truly-active TFs alone " +
                    "produced λ≈1.4-1.5 with a plain estimate); the truncated version narrows that to " +
                    "roughly λ≈1.0-1.1 for the same case, while still correctly recovering GENUINE " +
                    "inflation when it is really there. It narrows the ambiguity, it does not remove it. If " +
                    "you have independent reason to expect widespread real TF activity, read the PARAMETRIC " +
                    "column instead and treat the calibrated one as a conservative lower bound. See the README."));

                const tablePane = buildTable(block.table, {
                    maxHeight: 380,
                    rowClassFn: (row) => {
                        const fdr = row["FDR (calibrated)"];
                        return (fdr !== "" && fdr !== null && isFinite(Number(fdr)) && Number(fdr) < alpha) ? "row-significant" : null;
                    },
                    cellFormatFn: (col, val) => {
                        if (col === "p-value (parametric)" || col === "FDR (parametric)" ||
                            col === "p-value (permutation)" || col === "FDR (permutation)" ||
                            col === "p-value (calibrated)" || col === "FDR (calibrated)") return _fmtSci(val, 3);
                        if (col === "Activity Score" || col === "Calibrated Score (z)") return _fmtFixed(val, 3);
                        return val;
                    },
                });
                activityHost.appendChild(tablePane);
                const dlBtn = el("button", { class: "btn" }, "Download result table");
                dlBtn.addEventListener("click", () => {
                    downloadFile("/api/transcriptomics/regulon/download-activity?key=" + encodeURIComponent(db), "transcriptomics_regulon_activity_" + db.toLowerCase() + ".csv");
                });
                activityHost.appendChild(dlBtn);

                // QQ plot lives OUTSIDE the hasSig gate below: it explains
                // WHY a run has few (or no) significant TFs after
                // calibration, so it needs to render even when hasSig is
                // false - that is precisely the case it is most useful in.
                const qqControls = el("div", { class: "chart-toggle-row" });
                const qqArea = el("div", { class: "chart-pair" });
                activityHost.appendChild(qqControls);
                activityHost.appendChild(qqArea);
                const qqId = _oraChartId("calibrationqq", slug, null);
                const qqDiv = el("div", { id: qqId, class: "chart-container", style: { display: "none" } });
                qqArea.appendChild(qqDiv);
                let qqShown = false, qqBuilt = false;
                const qqBtn = el("button", { class: "btn chart-toggle-btn" }, "Show calibration QQ plot");
                qqBtn.addEventListener("click", () => {
                    if (!qqShown) {
                        qqDiv.style.display = "";
                        if (!qqBuilt) { buildRegulonCalibrationQQPlot(qqDiv, block.table, block, slug, {}); qqBuilt = true; }
                        else if (typeof Plotly !== "undefined") Plotly.Plots.resize(qqDiv);
                        qqShown = true; qqBtn.textContent = "Hide calibration QQ plot";
                    } else { qqDiv.style.display = "none"; qqShown = false; qqBtn.textContent = "Show calibration QQ plot"; }
                });
                qqControls.appendChild(qqBtn);

                const hasSig = (block.table.rows || []).some(row =>
                    row["FDR (calibrated)"] !== null && row["FDR (calibrated)"] !== "" &&
                    isFinite(Number(row["FDR (calibrated)"])) && Number(row["FDR (calibrated)"]) < alpha
                );
                if (!hasSig) return;

                const controls = el("div", { class: "chart-toggle-row" });
                const area = el("div", { class: "chart-pair" });
                activityHost.appendChild(controls);
                activityHost.appendChild(area);
                let topN = 25;
                let fontSize = 12;
                const barId = _oraChartId("activitybar", slug, null);
                const barDiv = el("div", { id: barId, class: "chart-container", style: { display: "none" } });
                const volcanoId = _oraChartId("activityvolcano", slug, null);
                const volcanoDiv = el("div", { id: volcanoId, class: "chart-container", style: { display: "none" } });
                const vsExprId = _oraChartId("activityvsexpr", slug, null);
                const vsExprDiv = el("div", { id: vsExprId, class: "chart-container", style: { display: "none" } });
                // On-screen key for the activity-vs-own-expression figure.
                // The READMEs are not reachable from the app yet, so a figure
                // whose meaning lives only there cannot be read at all - this
                // one uses four channels (fill, colour, ring, size) and was
                // demonstrably not self-explanatory without it.
                const vsExprLegend = el("div", { class: "caption chart-encoding-key", style: { display: "none" } }, [
                    el("div", {}, "How to read this figure:"),
                    el("div", {}, "• FILLED = the TF's ACTIVITY is significant (red = inferred active, blue = inferred repressed). Hollow grey = not significant."),
                    el("div", {}, "• DARK RING = the TF's OWN mRNA change is significant (its own padj passes the threshold). No ring = its own expression did not significantly change."),
                    el("div", {}, "• SIZE = how many of that TF's targets are present in your data (circle AREA tracks the count)."),
                    el("div", {}, "• A filled red circle with NO ring, sitting near the horizontal zero line, is the case worth noticing: the TF appears to be acting, but its own transcript never moved - the signature of activation after translation (phosphorylation, localisation, degradation) rather than by making more of the TF."),
                ]);
                area.appendChild(barDiv);
                area.appendChild(volcanoDiv);
                area.appendChild(vsExprLegend);
                area.appendChild(vsExprDiv);

                function redrawBar() { buildRegulonActivityBarPlot(barDiv, block.table, slug, alpha, { topN, fontSize }); }
                function redrawVolcano() { buildRegulonActivityVolcano(volcanoDiv, block.table, slug, alpha, { fontSize }); }
                function redrawVsExpr() { buildRegulonActivityVsExpression(vsExprDiv, block.table, geneLog2FC, genePadj, slug, alpha, { fontSize }); }

                let barShown = false, barBuilt = false;
                const barBtn = el("button", { class: "btn chart-toggle-btn" }, "Show TF activity bar plot");
                barBtn.addEventListener("click", () => {
                    if (!barShown) {
                        barDiv.style.display = "";
                        if (!barBuilt) { redrawBar(); barBuilt = true; } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(barDiv);
                        barShown = true; barBtn.textContent = "Hide TF activity bar plot";
                    } else { barDiv.style.display = "none"; barShown = false; barBtn.textContent = "Show TF activity bar plot"; }
                });
                controls.appendChild(barBtn);

                let volcanoShown = false, volcanoBuilt = false;
                const volcanoBtn = el("button", { class: "btn chart-toggle-btn" }, "Show TF activity volcano");
                volcanoBtn.addEventListener("click", () => {
                    if (!volcanoShown) {
                        volcanoDiv.style.display = "";
                        if (!volcanoBuilt) { redrawVolcano(); volcanoBuilt = true; } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(volcanoDiv);
                        volcanoShown = true; volcanoBtn.textContent = "Hide TF activity volcano";
                    } else { volcanoDiv.style.display = "none"; volcanoShown = false; volcanoBtn.textContent = "Show TF activity volcano"; }
                });
                controls.appendChild(volcanoBtn);

                let vsExprShown = false, vsExprBuilt = false;
                const vsExprBtn = el("button", { class: "btn chart-toggle-btn" }, "Show activity vs. own expression");
                vsExprBtn.addEventListener("click", () => {
                    if (!vsExprShown) {
                        vsExprDiv.style.display = "";
                        vsExprLegend.style.display = "";
                        if (!vsExprBuilt) { redrawVsExpr(); vsExprBuilt = true; } else if (typeof Plotly !== "undefined") Plotly.Plots.resize(vsExprDiv);
                        vsExprShown = true; vsExprBtn.textContent = "Hide activity vs. own expression";
                    } else {
                        vsExprDiv.style.display = "none";
                        vsExprLegend.style.display = "none";
                        vsExprShown = false; vsExprBtn.textContent = "Show activity vs. own expression";
                    }
                });
                controls.appendChild(vsExprBtn);

                const topNSel = el("label", { class: "chart-option-control" },
                    [el("span", { class: "chart-option-label" }, "TFs shown (bar plot)")]);
                const topNSelect = el("select", { class: "chart-option-select" });
                [["10", "10"], ["25", "25"], ["50", "50"], ["All significant", "all"]].forEach(([lbl, val]) =>
                    topNSelect.appendChild(el("option", { value: val }, lbl)));
                topNSelect.value = "25";
                topNSelect.addEventListener("change", () => {
                    topN = topNSelect.value === "all" ? Infinity : Number(topNSelect.value);
                    if (barShown) redrawBar();
                });
                topNSel.appendChild(topNSelect);
                controls.appendChild(topNSel);
            }
            dbSelect.addEventListener("change", () => renderActivityDb(dbSelect.value));
            renderActivityDb(activityDbs[0]);
        }
    }

    // --- Pane: Regulon Detail (per-TF) ---------------------------------------
    // Only TFs with a stored tf_regression_detail entry are selectable here -
    // see run_regulon_activity's docstring (REGULON_MAX_CURVES_STORED top TFs
    // by parametric FDR per database only, not every TF tested).
    if (activity) {
        const detailDbs = Object.keys(activity.tables).filter((db) => {
            const block = activity.tables[db];
            return block && block.tf_regression_detail && Object.keys(block.tf_regression_detail).length > 0;
        });
        if (detailDbs.length > 0) {
            const detailPane = nav.pane("Regulon Detail");
            detailPane.appendChild(el("div", { class: "caption" },
                "Shows the actual regression behind one TF's activity score: each of its targets plotted by " +
                "whether it's predicted to be activated or repressed, against how much that gene actually " +
                "moved in your data. The grey lines are what the SAME fit looks like on random relabelings of " +
                "which genes count as this TF's targets - a visual comparison, not a stronger proof than the " +
                "permutation p-value already in the table (both share the same known limitation for truly " +
                "co-regulated targets - see the caption on the TF Activity pane)."));

            const dbSelRow = el("div", { class: "row" });
            const detailDbSelect = el("select", { class: "text-input" });
            detailDbs.forEach((db) => detailDbSelect.appendChild(el("option", { value: db }, db)));
            dbSelRow.appendChild(el("div", { class: "col" }, detailDbSelect));
            const tfSelect = el("select", { class: "text-input" });
            dbSelRow.appendChild(el("div", { class: "col" }, tfSelect));
            detailPane.appendChild(dbSelRow);

            // One stable container div reused across TF selection changes -
            // mirrors GSEA's "Per gene set" pattern (renderTranscriptomicsGseaResultsScreen).
            const detailArea = el("div", { class: "chart-pair" });
            detailPane.appendChild(detailArea);
            const regressionDiv = el("div", { id: "chart-regulon-regression", class: "chart-container" });
            const coherenceDiv = el("div", { id: "chart-regulon-coherence", class: "chart-container" });
            detailArea.appendChild(regressionDiv);
            detailArea.appendChild(coherenceDiv);

            function renderDetailTf() {
                const db = detailDbSelect.value;
                const tf = tfSelect.value;
                const block = activity.tables[db];
                if (!block || !tf) return;
                const row = (block.table.rows || []).find((r) => r["TF"] === tf);
                const detail = (block.tf_regression_detail || {})[tf];
                if (!row || !detail) return;
                const nullSlopes = (block.null_slopes_by_signature || {})[detail.signature] || [];
                buildRegulonUlmRegressionPlot(regressionDiv, row, detail, nullSlopes, block, tf);
                buildRegulonCoherenceBars(coherenceDiv, row, detail, geneLog2FC, tf);
            }

            function populateTfSelect() {
                const db = detailDbSelect.value;
                const block = activity.tables[db];
                tfSelect.innerHTML = "";
                if (!block) return;
                Object.keys(block.tf_regression_detail || {}).forEach((tf) => {
                    tfSelect.appendChild(el("option", { value: tf }, tf));
                });
                renderDetailTf();
            }
            detailDbSelect.addEventListener("change", populateTfSelect);
            tfSelect.addEventListener("change", renderDetailTf);
            populateTfSelect();
        }
    }

    // --- Pane: TF Target Enrichment ------------------------------------------
    if (enrichment) {
        const enrichPane = nav.pane("TF Target Enrichment");
        const enrichKeys = Object.keys(enrichment.tables);
        if (enrichKeys.length === 0) {
            enrichPane.appendChild(alertBox("info", "Target enrichment was not run for this analysis."));
        } else {
            const selRow = el("div", { class: "row" });
            const keySelect = el("select", { class: "text-input" });
            enrichKeys.sort().forEach((key) => {
                const block = enrichment.tables[key];
                const suffix = block.n_terms_tested > 0 ? ` (${block.n_terms_tested} TFs)` : " (no TFs)";
                keySelect.appendChild(el("option", { value: key }, key + suffix));
            });
            selRow.appendChild(el("div", { class: "col" }, keySelect));
            enrichPane.appendChild(selRow);
            const enrichHost = el("div", {});
            enrichPane.appendChild(enrichHost);

            function renderEnrichKey(key) {
                enrichHost.innerHTML = "";
                const block = enrichment.tables[key];
                if (!block) return;
                const slug = _transcriptomicsRegulonSlug(key);
                if (block.note) enrichHost.appendChild(alertBox("info", block.note));
                if (typeof block.universe_size === "number" && block.n_terms_tested > 0) {
                    enrichHost.appendChild(metricCards([
                        ["TFs tested", block.n_terms_tested],
                        ["Universe size", block.universe_size],
                        ["DEGs in universe", block.n_deg_in_universe],
                    ], 3));
                }
                if (!block.table.rows || block.table.rows.length === 0) return;

                const tableEl = buildTable(block.table, {
                    maxHeight: 380,
                    rowClassFn: (row) => {
                        const fdr = row["FDR"];
                        return (fdr !== "" && fdr !== null && isFinite(Number(fdr)) && Number(fdr) < alpha) ? "row-significant" : null;
                    },
                    cellFormatFn: (col, val) => {
                        if (col === "p-value" || col === "FDR" || col === "Bonferroni") return _fmtSci(val, 3);
                        if (col === "Fold Enrichment") return _fmtFixed(val, 2);
                        return val;
                    },
                });
                enrichHost.appendChild(tableEl);
                const dlBtn = el("button", { class: "btn" }, "Download result table");
                dlBtn.addEventListener("click", () => {
                    // Per-key filename - the regulon enrichment screen has one
                    // table per database x direction, which all collided on a
                    // single filename before (same bug as ORA/GSEA above).
                    downloadFile("/api/transcriptomics/regulon/download?key=" + encodeURIComponent(key),
                                 "transcriptomics_regulon_enrichment_" + _transcriptomicsOraSlug(key) + ".csv");
                });
                enrichHost.appendChild(dlBtn);

                const hasSig = (block.table.rows || []).some(row =>
                    row["FDR"] !== null && row["FDR"] !== "" && isFinite(Number(row["FDR"])) && Number(row["FDR"]) < alpha
                );
                if (!hasSig) return;

                let topN = 20, sortBy = "fdr", fontSize = 12;
                let dotColor = ORA_DOT_COLOR_SCALES[0][1];
                let barColor = ORA_BAR_COLORS[0][1];

                const controls = el("div", { class: "chart-toggle-row" });
                const area = el("div", { class: "chart-pair" });
                enrichHost.appendChild(controls);
                enrichHost.appendChild(area);

                function addOpt(label, options, initial, onChange) {
                    const control = el("label", { class: "chart-option-control" });
                    control.appendChild(el("span", { class: "chart-option-label" }, label));
                    const sel = el("select", { class: "chart-option-select" });
                    options.forEach(([l, v]) => sel.appendChild(el("option", { value: v }, l)));
                    sel.value = initial;
                    sel.addEventListener("change", () => onChange(sel.value));
                    control.appendChild(sel);
                    controls.appendChild(control);
                    return sel;
                }
                function addToggle(kind, label, buildFn) {
                    const chartId = _oraChartId(kind, slug, null);
                    const chartDiv = el("div", { id: chartId, class: "chart-container", style: { display: "none" } });
                    const btn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
                    let shown = false, built = false;
                    function renderChart() { buildFn(chartDiv); built = true; }
                    btn.addEventListener("click", () => {
                        if (!shown) {
                            chartDiv.style.display = "";
                            if (!built) renderChart(); else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                            shown = true; btn.textContent = "Hide " + label;
                        } else { chartDiv.style.display = "none"; shown = false; btn.textContent = "Show " + label; }
                    });
                    controls.appendChild(btn);
                    area.appendChild(chartDiv);
                    return { redraw: () => { if (shown) renderChart(); } };
                }

                addOpt("TFs shown", [["10", "10"], ["20", "20"], ["50", "50"], ["All significant", "all"]], String(topN === Infinity ? "all" : topN), (v) => {
                    topN = v === "all" ? Infinity : Number(v);
                    dotToggle.redraw(); barToggle.redraw();
                });
                addOpt("Rank & show by", [["FDR", "fdr"], ["Fold enrichment", "fold"], ["Gene count", "k"], ["Gene ratio", "ratio"]], sortBy, (v) => {
                    sortBy = v;
                    dotToggle.redraw(); barToggle.redraw();
                });
                const dotToggle = addToggle("dot", "dot plot", (div) => buildDotPlot(div, block.table, slug, null, dotColor, { topN, sortBy, fontSize, alpha }));
                addOpt("Dot gradient", ORA_DOT_COLOR_SCALES.map(([l, v]) => [l, v]), dotColor, (v) => { dotColor = v; dotToggle.redraw(); });
                const barToggle = addToggle("bar", "bar plot", (div) => buildBarPlot(div, block.table, slug, null, barColor, { topN, sortBy, fontSize, alpha }));
                addOpt("Bar color", ORA_BAR_COLORS, barColor, (v) => { barColor = v; barToggle.redraw(); });

                // Figures C/D: regulons overlap heavily (e.g. RELA/NFKB1 can
                // share ~100 targets), so without this a user reports one
                // regulatory programme as several independent "findings".
                // Both reuse the SAME table shape ORA uses (RESULT_COLUMNS),
                // so the existing builders work completely unchanged.
                let edgeThreshold = 1;
                const enrichMapToggle = addToggle("enrichmap", "TF-TF overlap map", (div) =>
                    buildEnrichmentMap(div, block.table, slug, null, dotColor, { topN, sortBy, fontSize, edgeThreshold, alpha }));
                // "shared DE targets", not "shared targets": buildEnrichmentMap
                // draws edges from the Overlapping Genes column, which holds
                // each TF's targets that are DIFFERENTIALLY EXPRESSED here -
                // not its full recorded regulon.
                addOpt("Edge threshold (shared DE targets)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (v) => {
                    edgeThreshold = Number(v);
                    enrichMapToggle.redraw();
                });
                const networkToggle = addToggle("network", "TF-target network", (div) =>
                    buildTranscriptomicsGeneConceptNetwork(div, block.table, slug, geneLog2FC, { topN, sortBy, fontSize, alpha, dotColorScale: dotColor }));
            }
            keySelect.addEventListener("change", () => renderEnrichKey(keySelect.value));
            const sortedEnrichKeys = enrichKeys.sort();
            if (sortedEnrichKeys.length > 0) renderEnrichKey(sortedEnrichKeys[0]);
        }
    }

    nav.finish();

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.transcriptomics.regulonConfig = null;
        STATE.transcriptomics.regulonResults = null;
        STATE.transcriptomics.regulonRunLocked = false;
        STATE.transcriptomics._regulonDraft = null;
        STATE.aaPage = "transcriptomics_regulon_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// TF / Regulon chart builders (Activity bar plot, Activity volcano). Reuse
// _newOraPlot/_oraChartId/_magnitudeMarkerSizes/_wrapOraLabel/
// _oraCategoricalChartHeight so these charts export (SVG/PNG/TIFF) and behave
// identically to every other Omiscera chart, without redefining any of that
// plumbing.
// ---------------------------------------------------------------------------

function _regulonActivityRows(table, alpha) {
    return _tableToObjects(table).filter(r =>
        r["FDR (calibrated)"] !== null && r["FDR (calibrated)"] !== "" &&
        isFinite(Number(r["FDR (calibrated)"])) && Number(r["FDR (calibrated)"]) < alpha
    );
}

function buildRegulonActivityBarPlot(div, table, slug, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    let rows = _regulonActivityRows(table, alpha).sort((a, b) => Math.abs(Number(b["Activity Score"])) - Math.abs(Number(a["Activity Score"])));
    const topN = opts.topN || 25;
    if (isFinite(topN)) rows = rows.slice(0, topN);
    if (!rows.length) { div.textContent = "No significant TFs to plot."; return; }
    // Ascending by score so the diverging bars read bottom(repressed)-to-top(active).
    rows = rows.slice().sort((a, b) => Number(a["Activity Score"]) - Number(b["Activity Score"]));
    const names = rows.map(r => String(r["TF"] || ""));
    const scores = rows.map(r => Number(r["Activity Score"]));
    const chartHeight = _oraCategoricalChartHeight(names, fontSize);
    const trace = {
        x: scores,
        y: names,
        type: "bar",
        orientation: "h",
        marker: { color: scores.map(s => s >= 0 ? "#c0392b" : "#2874a6") },
        customdata: rows.map(r => [Number(r["Signed Targets Tested"]), Number(r["FDR (calibrated)"])]),
        hovertemplate: "<b>%{y}</b><br>Activity score: %{x:.3f}<br>Signed targets tested: %{customdata[0]}<br>FDR (calibrated): %{customdata[1]:.3g}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Activity score (red = active, blue = repressed)", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#999" },
        yaxis: { tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 160, r: 40, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "activitybar", "TF activity bar plot");
}

function buildRegulonActivityVolcano(div, table, slug, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : 0.05;
    const rows = _tableToObjects(table).filter(r => r["Activity Score"] !== null && isFinite(Number(r["Activity Score"])));
    if (!rows.length) { div.textContent = "No TFs to plot."; return; }
    const x = rows.map(r => Number(r["Activity Score"]));
    const y = rows.map(r => -Math.log10(Number(r["FDR (calibrated)"]) + 1e-300));
    const sizes = _magnitudeMarkerSizes(rows.map(r => Number(r["Signed Targets Tested"])), 8, 26);
    const colors = rows.map(r => Number(r["FDR (calibrated)"]) < sigAlpha ? (Number(r["Activity Score"]) >= 0 ? "#c0392b" : "#2874a6") : "#b8c2cc");
    const trace = {
        x, y,
        mode: "markers",
        type: "scatter",
        text: rows.map(r => r["TF"]),
        customdata: rows.map(r => [Number(r["Signed Targets Tested"])]),
        marker: { color: colors, size: sizes, opacity: 0.85 },
        hovertemplate: "<b>%{text}</b><br>Activity score: %{x:.3f}<br>-log10(FDR, calibrated): %{y:.2f}<br>Signed targets tested: %{customdata[0]}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Activity score", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { title: { text: "-log10(FDR, calibrated)", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: false },
        height: 500,
        margin: { l: 95, r: 40, t: 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "activityvolcano", "TF activity volcano");
}


// Inverse standard normal CDF (probit), Acklam's rational approximation -
// relative error < 1.15e-9, plenty for a visual QQ plot. No JS stdlib
// equivalent exists, and the calibration model itself treats scores via a
// normal approximation (ndtr, server-side) - see _empirical_null_calibration's
// docstring - so a normal-quantile QQ plot is the correct diagnostic to pair
// with it, not an approximation error introduced here.
function _probit(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
               6.680131188771972e+01, -1.328068155288572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
               3.754408661907416e+00];
    const pLow = 0.02425, pHigh = 1 - pLow;
    let q, r;
    if (p < pLow) {
        q = Math.sqrt(-2 * Math.log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p <= pHigh) {
        q = p - 0.5;
        r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    } else {
        q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
}

// ---------------------------------------------------------------------------
// Calibration QQ plot — the figure that SHOWS lambda is real rather than
// asking the user to trust a number. Observed TF activity scores (sorted)
// against the quantiles a perfectly-calibrated null (mean 0, sd 1) would
// predict. A run with no inflation problem sits on y=x; the fitted empirical-
// null line (intercept = null_center, slope = lambda_used) shows what this
// run's ACTUAL null looks like - visibly steeper than y=x whenever lambda > 1.
// ---------------------------------------------------------------------------
function buildRegulonCalibrationQQPlot(div, table, block, slug, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const rows = _tableToObjects(table).filter(r => r["Activity Score"] !== null && isFinite(Number(r["Activity Score"])));
    if (rows.length < 5) { div.textContent = "Too few TFs tested for a calibration QQ plot."; return; }
    const scores = rows.map(r => Number(r["Activity Score"])).sort((a, b) => a - b);
    const nPts = scores.length;
    const theoretical = scores.map((_, i) => _probit((i + 0.5) / nPts));

    // Manual min/max (not Math.min(...arr)): nPts can run into the
    // thousands, and spreading an array that large into Math.min/max risks a
    // stack overflow - same rationale _magnitudeMarkerSizes already applies.
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < nPts; i++) {
        if (theoretical[i] < lo) lo = theoretical[i];
        if (theoretical[i] > hi) hi = theoretical[i];
        if (scores[i] < lo) lo = scores[i];
        if (scores[i] > hi) hi = scores[i];
    }
    const pad = (hi - lo) * 0.05 || 1;
    lo -= pad; hi += pad;

    const lambdaUsed = Number(block.lambda_used) || 1.0;
    const nullCenter = Number(block.null_center) || 0.0;

    const traces = [
        {
            x: theoretical, y: scores, mode: "markers", type: "scatter", name: "Observed TF scores",
            marker: { color: "#1B639C", size: 5, opacity: 0.6 },
            hovertemplate: "theoretical: %{x:.2f}<br>observed: %{y:.2f}<extra></extra>",
        },
        {
            x: [lo, hi], y: [lo, hi], mode: "lines", type: "scatter", name: "y = x (no inflation)",
            line: { color: "#999", width: 1, dash: "dot" }, hoverinfo: "skip",
        },
        {
            x: [lo, hi], y: [nullCenter + lambdaUsed * lo, nullCenter + lambdaUsed * hi],
            mode: "lines", type: "scatter", name: "Fitted empirical null (slope = λ)",
            line: { color: "#c0392b", width: 2 }, hoverinfo: "skip",
        },
    ];
    const layout = {
        font: { size: fontSize },
        title: { text: "λ = " + lambdaUsed.toFixed(3), font: { size: 14 } },
        xaxis: { title: { text: "Theoretical quantile (no inflation)", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, range: [lo, hi] },
        yaxis: { title: { text: "Observed activity score", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, range: [lo, hi] },
        height: 480, margin: { l: 90, r: 30, t: 46, b: 56 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        legend: { orientation: "h", y: -0.15 },
    };
    _newOraPlot(div, traces, layout, slug, null, "calibrationqq", "Calibration QQ plot");
}


// ---------------------------------------------------------------------------
// Figure B — Activity vs. the TF's own expression. Three numbers that already
// exist (activity score, the TF's own log2FC, the TF's own padj) - no new
// statistics. Quadrants are biologically real: active + own-mRNA-up =
// transcriptionally driven; active + own-mRNA-flat = post-translationally
// activated (phosphorylation, localization, degradation - the TF's mRNA level
// doesn't have to move for its protein activity to change).
//
// ENCODING, and why it changed. An earlier version used marker FILL for the
// activity FDR and nothing at all for the TF's own padj. Because a large
// |activity score| is exactly what produces a small activity FDR, fill was
// almost perfectly redundant with x-position - the strongest visual channel
// in the figure repeated the x-axis while the y-axis carried magnitude with
// no confidence attached. Worse, it made the figure unable to show the very
// case it exists for: the README promised "strongly active with an open
// marker = post-translational", but under that encoding "strongly active"
// implies a significant activity FDR, hence FILLED, so that reading could
// essentially never occur. Now:
//
//   FILL / COLOUR : the ACTIVITY call is significant (calibrated FDR < alpha);
//                   red = inferred active, blue = inferred repressed.
//   DARK RING     : the TF's OWN mRNA change is significant (own padj < alpha).
//   SIZE          : signed targets tested (area-proportional, via _magnitudeMarkerSizes).
//
// So a filled red marker with NO ring sitting near y = 0 is the money case:
// the TF is acting, but its own transcript never moved.
// ---------------------------------------------------------------------------
function buildRegulonActivityVsExpression(div, table, geneLog2FC, genePadj, slug, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    genePadj = genePadj || {};
    const fontSize = opts.fontSize || 12;
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : 0.05;
    const allRows = _tableToObjects(table).filter(r => r["Activity Score"] !== null && isFinite(Number(r["Activity Score"])));
    const rows = allRows.filter(r => Object.prototype.hasOwnProperty.call(geneLog2FC, String(r["TF"] || "").toUpperCase()));
    if (!rows.length) { div.textContent = "None of the tested TFs' own genes are present in the resolved DE table."; return; }
    const nOmitted = allRows.length - rows.length;
    const sym = (r) => String(r["TF"] || "").toUpperCase();
    const x = rows.map(r => Number(r["Activity Score"]));
    const y = rows.map(r => geneLog2FC[sym(r)]);
    const activitySig = rows.map(r => Number(r["FDR (calibrated)"]) < sigAlpha);
    // A TF whose own padj is missing is NOT treated as significant - absence
    // of evidence must not render as a confidence ring.
    const ownSig = rows.map((r) => {
        const p = genePadj[sym(r)];
        return p !== undefined && isFinite(p) && p < sigAlpha;
    });
    const sizes = _magnitudeMarkerSizes(rows.map(r => Number(r["Signed Targets Tested"])), 8, 26);
    const fills = rows.map((r, i) => activitySig[i] ? (x[i] >= 0 ? "#c0392b" : "#2874a6") : "#b8c2cc");
    const lineColors = ownSig.map((s) => s ? "#1a1a1a" : "rgba(0,0,0,0)");
    const lineWidths = ownSig.map((s) => s ? 2.2 : 0);
    const trace = {
        x, y,
        mode: "markers",
        type: "scatter",
        text: rows.map(r => r["TF"]),
        customdata: rows.map((r, i) => [
            activitySig[i] ? "significant" : "not significant",
            genePadj[sym(r)] === undefined ? "n/a" : Number(genePadj[sym(r)]).toExponential(2),
            ownSig[i] ? "significant" : "not significant",
            Number(r["Signed Targets Tested"]),
        ]),
        marker: {
            color: fills, size: sizes, opacity: 0.85, symbol: "circle",
            line: { width: lineWidths, color: lineColors },
        },
        hovertemplate:
            "<b>%{text}</b><br>Activity score: %{x:.3f} (%{customdata[0]})" +
            "<br>Own log2FC: %{y:.3f}<br>Own padj: %{customdata[1]} (%{customdata[2]})" +
            "<br>Signed targets tested: %{customdata[3]}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Activity score", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#999" },
        yaxis: { title: { text: "TF's own log2 fold-change", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#999" },
        height: 520,
        margin: { l: 95, r: 40, t: nOmitted > 0 ? 46 : 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
        annotations: nOmitted > 0 ? [{
            xref: "paper", yref: "paper", x: 0, y: 1.06, showarrow: false, xanchor: "left",
            text: nOmitted + " TF(s) not in the resolved DE table are omitted from this figure.",
            font: { size: Math.max(9, fontSize - 2), color: "#888" },
        }] : [],
    };
    _newOraPlot(div, [trace], layout, slug, null, "activityvsexpr", "Activity vs. own expression");
}


// ---------------------------------------------------------------------------
// Figure A (Regulon Detail pane) — the actual ULM regression, drawn. x =
// target weight (+1 activated / -1 repressed), y = the gene's ranking metric.
// The grey lines are the SAME fit recomputed on random relabelings (the null
// this module's permutation p-value is built from) - shown as a visual
// comparison, explicitly NOT framed as stronger evidence than the p-value
// already in the table (see this session's finding: gene-label permutation
// shares the parametric test's own blind spot for truly co-regulated
// targets, so the fan is illustrative, not a superior proof).
// ---------------------------------------------------------------------------
function buildRegulonUlmRegressionPlot(div, row, detail, nullSlopes, block, tf) {
    if (typeof Plotly === "undefined") return;
    // y-values come from the SERVER (detail.activated_values /
    // repressed_values) - these are the ranking-metric values the fit
    // actually consumed. They are NOT the client's log2FC lookup: the two
    // are different variables whenever the ranking metric isn't "log2fc"
    // (the default is signed -log10(p)), and an earlier version of this
    // figure plotted log2FC against a slope fitted in metric units, so the
    // line did not pass through its own points.
    const activated = detail.activated_symbols || [];
    const repressed = detail.repressed_symbols || [];
    const actY = detail.activated_values || [];
    const repY = detail.repressed_values || [];
    const metricLabel = GSEA_RANKING_METRIC_LABELS[block.metric] || block.metric || "Ranking metric";

    const traces = [];

    // x=0 mass: the non-target genes the fit is anchored by. Drawn as a
    // violin (not ~3000 markers) - see run_regulon_activity's
    // ranking_values_sample comment for why this is the all-genes
    // distribution and labelled as such.
    const allSample = block.ranking_values_sample || [];
    if (allSample.length) {
        traces.push({
            x: allSample.map(() => 0), y: allSample,
            type: "violin", name: "All tested genes (weight 0)",
            line: { color: "#9DB2C6" }, fillcolor: "rgba(157,178,198,0.35)",
            points: false, width: 0.9, hoverinfo: "skip", showlegend: true,
        });
    }

    traces.push({
        x: actY.map(() => 1 + (Math.random() - 0.5) * 0.3),
        y: actY, mode: "markers", type: "scatter", name: "Activated targets",
        marker: { color: "#c0392b", size: 7, opacity: 0.75 },
        text: activated, hovertemplate: "%{text}<br>" + metricLabel + ": %{y:.3f}<extra></extra>",
    });
    traces.push({
        x: repY.map(() => -1 + (Math.random() - 0.5) * 0.3),
        y: repY, mode: "markers", type: "scatter", name: "Repressed targets",
        marker: { color: "#2874a6", size: 7, opacity: 0.75 },
        text: repressed, hovertemplate: "%{text}<br>" + metricLabel + ": %{y:.3f}<extra></extra>",
    });

    // Every OLS line passes through the centroid (x_mean, y_mean), so only
    // the slope varies between the observed fit and each null draw. BOTH
    // come from the server's own fit - x_mean is (n_pos-n_neg)/n over ALL
    // tested genes (not over the targets), y_mean is the mean ranking value
    // over all genes.
    const xBar = Number(detail.x_mean) || 0;
    const yBar = Number(block.y_mean) || 0;
    const x0 = -1.3, x1 = 1.3;
    const lineY = (s) => [yBar + s * (x0 - xBar), yBar + s * (x1 - xBar)];

    const sampleCount = Math.min(100, nullSlopes.length);
    const step = nullSlopes.length > sampleCount ? Math.floor(nullSlopes.length / sampleCount) : 1;
    for (let i = 0; i < nullSlopes.length; i += step) {
        traces.push({
            x: [x0, x1], y: lineY(nullSlopes[i]),
            mode: "lines", type: "scatter", showlegend: false, hoverinfo: "skip",
            line: { color: "rgba(120,120,120,0.12)", width: 1 },
        });
    }

    const beta1 = Number(detail.beta1);
    traces.push({
        x: [x0, x1], y: lineY(beta1),
        mode: "lines", type: "scatter", name: "Fitted slope",
        line: { color: beta1 >= 0 ? "#c0392b" : "#2874a6", width: 3 },
        hoverinfo: "skip",
    });

    const layout = {
        title: { text: _escapePlotlyText(tf), font: { size: 15 } },
        xaxis: { title: { text: "Target weight (repressed ← → activated)", font: { size: 13 } },
                 tickvals: [-1, 0, 1], ticktext: ["Repressed", "Not a target", "Activated"], range: [-1.6, 1.6] },
        yaxis: { title: { text: metricLabel, font: { size: 13 } } },
        height: 480, margin: { l: 80, r: 30, t: 46, b: 56 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        legend: { orientation: "h", y: -0.15 },
        annotations: [{
            xref: "paper", yref: "paper", x: 0.02, y: 0.98, showarrow: false, xanchor: "left", yanchor: "top",
            text: "Activity score " + Number(row["Activity Score"]).toFixed(2) +
                  " | p (calibrated) " + Number(row["p-value (calibrated)"]).toExponential(2) +
                  " | p (parametric) " + Number(row["p-value (parametric)"]).toExponential(2) +
                  " | p (permutation) " + Number(row["p-value (permutation)"]).toExponential(2),
            font: { size: 11, color: "#555" },
            bgcolor: "rgba(255,255,255,0.7)",
        }],
    };
    _newOraPlot(div, traces, layout, "regulon_detail", tf, "ulmregression", "ULM regression");
}


// ---------------------------------------------------------------------------
// Figure E (Regulon Detail pane) — coherence: what fraction of this TF's
// activated targets actually moved UP, and what fraction of its repressed
// targets actually moved DOWN. Separates COHERENCE from MAGNITUDE - a TF can
// have a large activity score driven by a few big movers, or a modest score
// with near-unanimous agreement across many targets; the table alone can't
// show which.
// ---------------------------------------------------------------------------
function buildRegulonCoherenceBars(div, row, detail, geneLog2FC, tf) {
    if (typeof Plotly === "undefined") return;
    const activated = (detail.activated_symbols || []).map((s) => geneLog2FC[s]).filter((v) => v !== undefined && isFinite(v));
    const repressed = (detail.repressed_symbols || []).map((s) => geneLog2FC[s]).filter((v) => v !== undefined && isFinite(v));
    if (!activated.length && !repressed.length) { div.textContent = "No target log2FC values available."; return; }
    // >= 0 / < 0 (not > 0 / <= 0), to match build_ranked_list's own
    // sign(0) := +1 convention (Transcriptomics_module.py) - a gene with
    // log2FC exactly 0 is "up" to the ranking metric the activity score was
    // fitted on, so it must be "up" here too, or this figure could show a
    // coherence pattern that contradicts the table it explains. No such gene
    // exists in the bundled example, but the two must agree by construction,
    // not by the coincidence of a particular dataset.
    const actUp = activated.filter((v) => v >= 0).length;
    const actDown = activated.length - actUp;
    const repDown = repressed.filter((v) => v < 0).length;
    const repUp = repressed.length - repDown;

    const traces = [
        {
            y: ["Activated targets", "Repressed targets"],
            x: [activated.length ? actUp / activated.length : 0, repressed.length ? repUp / repressed.length : 0],
            name: "Moved UP", orientation: "h", type: "bar", marker: { color: "#c0392b" },
            customdata: [actUp, repUp], hovertemplate: "%{customdata} gene(s) moved up<extra></extra>",
        },
        {
            y: ["Activated targets", "Repressed targets"],
            x: [activated.length ? actDown / activated.length : 0, repressed.length ? repDown / repressed.length : 0],
            name: "Moved DOWN", orientation: "h", type: "bar", marker: { color: "#2874a6" },
            customdata: [actDown, repDown], hovertemplate: "%{customdata} gene(s) moved down<extra></extra>",
        },
    ];
    const layout = {
        title: { text: "Coherence: " + _escapePlotlyText(tf), font: { size: 14 } },
        barmode: "stack",
        xaxis: { title: { text: "Fraction of targets", font: { size: 12 } }, range: [0, 1], tickformat: ".0%" },
        yaxis: { automargin: true },
        height: 220, margin: { l: 130, r: 30, t: 40, b: 46 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        legend: { orientation: "h", y: -0.3 },
    };
    _newOraPlot(div, traces, layout, "regulon_detail", tf, "coherence", "Coherence bars");
}



// ---------------------------------------------------------------------------
// Transcriptomics — Network / Module Detection (fourth analysis tile)
// Seeds the STRING interaction network with the DE genes, propagates, then
// clusters the significant hits into modules. Reuses the Genomics netprop
// figures (buildNetpropNetwork / buildNetpropHistogram) unchanged - the server
// payload is deliberately the same shape - and the ORA dot/bar builders for
// per-module enrichment via the shared RESULT_COLUMNS contract.
// ---------------------------------------------------------------------------

// Mirrors Transcriptomics_network.NETWORK_N_PERMUTATIONS_CHOICES - 200 is
// deliberately NOT offered (see that constant's own comment: measured to
// underestimate the null variance badly enough to report false extreme
// significance on real data).
const TX_NETWORK_PERMUTATION_CHOICES = [["500", "500"], ["1000", "1000"], ["2000", "2000"]];
const TX_NETWORK_DIRECTIONS = [["All DE genes", "all"], ["Up-regulated only", "up"], ["Down-regulated only", "down"]];
const TX_NETWORK_WEIGHTINGS = [
    ["Equal weight (recommended)", "equal"],
    ["Weighted by |ranking metric|", "abs_metric"],
];
const TX_NETWORK_ENRICHMENT_DBS = [["GO", "GO"], ["MSigDB", "MSigDB"]];

// Rough wall-clock guide so a user isn't surprised by a multi-minute run.
// Anchored on measurements across the full offered range (500/1000/2000
// permutations over the human network took ~154/301/582 s - see
// tests/calibration_network_production_scale.py, Part B): cost is very
// close to linear in the permutation count, ~0.3 s/permutation.
function _txNetworkRuntimeHint(nPermutations) {
    const seconds = Math.round((Number(nPermutations) || 1000) * 0.3);
    if (seconds < 90) return "about " + seconds + " seconds";
    return "about " + Math.round(seconds / 60) + " minute(s)";
}

function renderTranscriptomicsNetworkConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" },
        "Transcriptomics — Network / Module Detection: setup"));

    const locked = STATE.transcriptomics.networkRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    if (!STATE.transcriptomics._networkDraft) {
        STATE.transcriptomics._networkDraft = {
            confidence: 700,
            direction: "all",
            weighting: "equal",
            metric: "signed_log_pvalue",
            padjCutoff: 0.05,
            lfcCutoff: 1.0,
            fdrAlpha: 0.05,
            nPermutations: 1000,
            moduleDetection: true,
            minModuleSize: 3,
            maxHitsForModules: 600,
            enrichmentDatabase: "GO",
            msigdbCollections: [],
            includeIea: true,
        };
    }
    const draft = STATE.transcriptomics._networkDraft;
    const taxonId = STATE.transcriptomics.taxonId;

    wrap.appendChild(el("div", { class: "caption" },
        "Organism: " + (STATE.transcriptomics.taxonName || "?") + " (STRING taxon " + (taxonId || "?") + "). " +
        "This analysis spreads your differential-expression signal through the STRING protein-interaction " +
        "network, so it can surface genes you did NOT upload but that sit at the centre of the response."));

    // --- seeding ---
    wrap.appendChild(el("div", { class: "section-label" }, "Seed genes"));
    const seedBox = el("div", { class: "rank-card ora-db-card" });

    function selectRow(box, label, choices, current, onChange) {
        const row = el("div", { class: "ora-term-size-row" });
        row.appendChild(el("span", { class: "chart-option-label" }, label));
        const sel = el("select", { class: "chart-option-select" });
        choices.forEach(([lbl, val]) => sel.appendChild(el("option", { value: val }, lbl)));
        sel.value = String(current);
        sel.disabled = locked;
        sel.addEventListener("change", () => onChange(sel.value));
        row.appendChild(sel);
        box.appendChild(row);
        return sel;
    }
    function numberRow(box, label, value, step, min, onChange) {
        const row = el("div", { class: "ora-term-size-row" });
        row.appendChild(el("span", { class: "chart-option-label" }, label));
        const inp = el("input", { type: "number", class: "ora-term-size-input", step: String(step), min: String(min) });
        inp.value = value;
        inp.disabled = locked;
        inp.addEventListener("change", () => onChange(inp.value));
        row.appendChild(inp);
        box.appendChild(row);
        return inp;
    }

    numberRow(seedBox, "Adjusted p-value cutoff (padj <=)", draft.padjCutoff, 0.001, 0, (v) => {
        const x = parseFloat(v); draft.padjCutoff = (isFinite(x) && x > 0 && x <= 1) ? x : 0.05;
    });
    numberRow(seedBox, "Minimum |log2 fold-change|", draft.lfcCutoff, 0.1, 0, (v) => {
        const x = parseFloat(v); draft.lfcCutoff = (isFinite(x) && x >= 0) ? x : 1.0;
    });
    selectRow(seedBox, "Direction", TX_NETWORK_DIRECTIONS, draft.direction, (v) => { draft.direction = v; });
    const metricRow = el("div", {});
    selectRow(seedBox, "Seed weighting", TX_NETWORK_WEIGHTINGS, draft.weighting, (v) => {
        draft.weighting = v;
        metricRow.classList.toggle("hidden", v !== "abs_metric");
    });
    selectRow(metricRow, "Ranking metric", Object.keys(GSEA_RANKING_METRIC_LABELS).map(
        (k) => [GSEA_RANKING_METRIC_LABELS[k], k]), draft.metric, (v) => { draft.metric = v; });
    metricRow.classList.toggle("hidden", draft.weighting !== "abs_metric");
    seedBox.appendChild(metricRow);
    seedBox.appendChild(el("div", { class: "caption ora-custom-help" },
        "Only genes passing BOTH cutoffs seed the network. \"Equal weight\" is recommended: the cutoffs have " +
        "already made the evidence call, so re-weighting by the same p-values counts it twice. Direction is " +
        "applied by choosing WHICH genes seed - never as a positive/negative weight, because the random walk " +
        "spreads an unsigned quantity and a signed weight would silently bias the whole run."));
    wrap.appendChild(seedBox);

    // --- network + significance ---
    wrap.appendChild(el("div", { class: "section-label" }, "Network and significance"));
    const netBox = el("div", { class: "rank-card ora-db-card" });
    // Read at call time, not module load: this const is declared later in the file.
    selectRow(netBox, "STRING confidence cutoff", NETPROP_CONFIDENCE_CHOICES, draft.confidence,
        (v) => { draft.confidence = Number(v); });
    numberRow(netBox, "Significance threshold (FDR <)", draft.fdrAlpha, 0.001, 0, (v) => {
        const x = parseFloat(v); draft.fdrAlpha = (isFinite(x) && x > 0 && x <= 1) ? x : 0.05;
    });
    const runtimeCaption = el("div", { class: "caption ora-custom-help" },
        "Estimated run time: " + _txNetworkRuntimeHint(draft.nPermutations) + ".");
    selectRow(netBox, "Number of permutations", TX_NETWORK_PERMUTATION_CHOICES, draft.nPermutations, (v) => {
        draft.nPermutations = Number(v);
        runtimeCaption.textContent = "Estimated run time: " + _txNetworkRuntimeHint(draft.nPermutations) + ".";
    });
    netBox.appendChild(runtimeCaption);
    netBox.appendChild(el("div", { class: "caption ora-custom-help" },
        "STRING's combined score folds in text-mining and cross-organism transfer as well as experiments, so " +
        "it is not purely experimental evidence and leans toward well-studied genes - which is why 0.7 (high) " +
        "is the default. Each permutation re-runs the whole propagation, so the count drives the run time."));
    wrap.appendChild(netBox);

    // --- modules ---
    wrap.appendChild(el("div", { class: "section-label" }, "Modules"));
    const modBox = el("div", { class: "rank-card ora-db-card" });
    const modRow = el("div", { class: "ora-iea-row" });
    const modCb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
    modCb.checked = !!draft.moduleDetection;
    modCb.disabled = locked;
    const modOptions = el("div", {});
    modCb.addEventListener("change", () => {
        draft.moduleDetection = modCb.checked;
        modOptions.classList.toggle("hidden", !modCb.checked);
    });
    modRow.appendChild(modCb);
    modRow.appendChild(document.createTextNode(" Detect modules and label each by functional enrichment"));
    modBox.appendChild(modRow);

    numberRow(modOptions, "Minimum module size", draft.minModuleSize, 1, 2, (v) => {
        const x = parseInt(v, 10); draft.minModuleSize = (isFinite(x) && x >= 2) ? x : 3;
    });
    numberRow(modOptions, "Maximum hits used for clustering", draft.maxHitsForModules, 50, 50, (v) => {
        const x = parseInt(v, 10); draft.maxHitsForModules = (isFinite(x) && x >= 50) ? x : 600;
    });
    const msigRow = el("div", {});
    selectRow(modOptions, "Label modules using", TX_NETWORK_ENRICHMENT_DBS, draft.enrichmentDatabase, (v) => {
        draft.enrichmentDatabase = v;
        msigRow.classList.toggle("hidden", v !== "MSigDB");
        renderRunBtn();
    });
    // Same helper (and same draft.msigdbCollections field name) ORA and GSEA
    // use, so the collection list and its per-organism gating stay identical
    // across all three without a second implementation.
    msigRow.appendChild(renderMsigdbCollectionSelector(locked, draft, taxonId, () => renderRunBtn()));
    msigRow.classList.toggle("hidden", draft.enrichmentDatabase !== "MSigDB");
    modOptions.appendChild(msigRow);
    modOptions.appendChild(el("div", { class: "caption ora-custom-help" },
        "Modules are labelled against the whole network's genes, not just the genes you uploaded - otherwise " +
        "every recruited gene would silently vanish from its own module's label. Clustering uses the strongest " +
        "hits by Z-score up to the cap, which keeps the grouping stable: on a very large hit set the algorithm " +
        "can otherwise collapse everything into one giant, uninformative \"module\". " +
        "Reactome is not offered here because its identifier bridge is built from your uploaded table, so " +
        "recruited genes would be dropped entirely."));
    modOptions.classList.toggle("hidden", !draft.moduleDetection);
    modBox.appendChild(modOptions);
    wrap.appendChild(modBox);

    // --- actions ---
    const actions = el("div", { class: "right-actions ora-config-actions" });
    const btnRow = el("div", { class: "btn-row ora-config-btn-row" });
    const backBtn = el("button", { class: "btn ora-back-btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "transcriptomics_module"; render(); });
    btnRow.appendChild(backBtn);
    const runWrap = el("div", { class: "col" });
    btnRow.appendChild(runWrap);
    actions.appendChild(btnRow);

    function renderRunBtn() {
        runWrap.innerHTML = "";
        const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
        const needsCollection = draft.enrichmentDatabase === "MSigDB" && draft.moduleDetection
            && (draft.msigdbCollections || []).length === 0;
        runBtn.disabled = locked || !taxonId || needsCollection;
        runBtn.addEventListener("click", () => startTranscriptomicsNetworkRun(draft));
        runWrap.appendChild(runBtn);
    }
    renderRunBtn();
    wrap.appendChild(actions);

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}

function startTranscriptomicsNetworkRun(draft) {
    const payload = {
        confidence: draft.confidence,
        direction: draft.direction,
        weighting: draft.weighting,
        metric: draft.metric,
        padj_cutoff: draft.padjCutoff,
        lfc_cutoff: draft.lfcCutoff,
        fdr_alpha: draft.fdrAlpha,
        n_permutations: draft.nPermutations,
        module_detection: !!draft.moduleDetection,
        min_module_size: draft.minModuleSize,
        max_hits_for_modules: draft.maxHitsForModules,
        enrichment_database: draft.enrichmentDatabase,
        msigdb_collections: draft.msigdbCollections || [],
        include_iea: !!draft.includeIea,
    };
    apiPostJson("/api/transcriptomics/network/configure", payload).then(() => {
        STATE.transcriptomics.networkConfig = payload;
        STATE.transcriptomics.networkRunLocked = true;
        STATE.transcriptomics.networkResults = null;
        // Navigate BEFORE starting the run so the progress bar divs exist by
        // the time pollJob's first tick arrives.
        STATE.aaPage = "transcriptomics_network_results";
        render();
        return apiPostJson("/api/transcriptomics/network/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("transcriptomics-network-progress-inner");
                const label = document.getElementById("transcriptomics-network-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => { STATE.transcriptomics.networkResults = result; render(); },
            onError: (err) => {
                alert("Network / Module Detection run failed: " + err);
                STATE.transcriptomics.networkRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.transcriptomics.networkRunLocked = false;
        render();
    });
}

function renderTranscriptomicsNetworkResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" },
        "Transcriptomics — Network / Module Detection: results"));

    const payload = STATE.transcriptomics.networkResults;
    if (!payload) {
        wrap.appendChild(alertBox("info", "Running analysis... (downloading the STRING network and running the permutation test can take a few minutes)"));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "transcriptomics-network-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "transcriptomics-network-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const prop = payload.propagation || {};
    const mapping = payload.mapping || {};
    const seeds = payload.seeds || {};
    const params = payload.parameters || {};
    const DB = "transcriptomics_network";   // keeps chart ids/filenames distinct from the Genomics tile

    (payload.caveats || []).forEach((c) => wrap.appendChild(alertBox("warning", c)));

    const card = el("div", { class: "rank-card ora-db-card ora-result-card" });
    wrap.appendChild(card);
    const navHost = el("div", {});
    card.appendChild(navHost);
    const nav = makeResultNav(navHost, { stateKey: "transcriptomics_network" });

    // ----- Overview -----
    const overview = nav.pane("Overview");
    overview.appendChild(metricCards([
        ["Genes seeded", seeds.n_seeds || 0],
        ["Significant", prop.n_hits || 0],
        ["Recruited", prop.n_recruited || 0],
        ["Modules", (payload.modules && (payload.modules.table.rows || []).length) || 0],
        ["Network proteins", prop.n_proteins || 0],
    ], 5));
    overview.appendChild(el("div", { class: "caption" },
        "\"Recruited\" genes are proteins you did NOT upload that the network places at the centre of your " +
        "response - the information this analysis adds over Over-Representation Analysis and GSEA. A seed " +
        "scoring significantly is largely expected: it receives signal directly and is compared against " +
        "permutations in which it is usually not a seed, so read the recruited set as the discovery."));
    const nd = prop.null_diagnostics || {};
    overview.appendChild(el("div", { class: "caption" },
        "Significance comes from " + (nd.n_permutations || params.n_permutations || "?") + " degree-preserving " +
        "permutations. p-values use a skew-aware (Pearson III) fit where the null is skewed enough to fit " +
        "reliably (" + (nd.n_gamma_tail || 0) + " of " + (prop.n_proteins || "?") + " proteins here), falling " +
        "back to the normal approximation otherwise. Measured null skewness at the top protein was " +
        (nd.null_skew_at_max !== undefined && nd.null_skew_at_max !== null ? Number(nd.null_skew_at_max).toFixed(2) : "n/a") +
        " (0 = perfectly normal); treat the ranking as more reliable than any single protein's exact FDR."));

    // ----- Gene mapping -----
    const mapPane = nav.pane("Gene mapping");
    mapPane.appendChild(metricCards([
        ["Uploaded genes", mapping.n_input || 0],
        ["Matched to STRING", mapping.n_mapped || 0],
        ["Not matched", mapping.n_unmapped || 0],
        ["Coverage", ((mapping.fraction_mapped || 0) * 100).toFixed(1) + "%"],
    ], 4));
    const tiers = mapping.n_by_tier || {};
    mapPane.appendChild(el("div", { class: "caption" },
        "Matched by gene symbol: " + (tiers.gene_symbol || 0) + ", by Ensembl ID: " + (tiers.ensembl_id || 0) +
        ", by Entrez ID: " + (tiers.entrez_id || 0) + ". " +
        (mapping.n_collisions_used ? mapping.n_collisions_used + " identifier(s) matched more than one protein; " +
            "the best-connected protein was used. " : "") +
        (seeds.n_multi_row_proteins ? seeds.n_multi_row_proteins + " seed protein(s) were mapped by MORE THAN " +
            "ONE row of your table (e.g. duplicate probes/transcripts for one gene); each still counts as a " +
            "single seed, weighted by its strongest evidence, not summed across rows. " : "") +
        "STRING names proteins using each organism's own nomenclature, so coverage varies by species."));
    if (mapping.table) mapPane.appendChild(buildTable(mapping.table, { maxHeight: 420 }));
    const mapDl = el("button", { class: "btn" }, "Download gene mapping (CSV)");
    mapDl.addEventListener("click", () => downloadFile(
        "/api/transcriptomics/network/download?which=mapping", "transcriptomics_network_gene_mapping.csv"));
    mapPane.appendChild(mapDl);

    // ----- Hit table -----
    const hitPane = nav.pane("Hit table");
    hitPane.appendChild(el("div", { class: "caption" },
        "Every protein in the network, ranked by FDR. \"Seed\" marks genes you uploaded that passed the " +
        "cutoffs; rows with no log2FC/padj were never in your table at all - those are the recruited ones."));
    const hitAlpha = params.fdr_alpha || 0.05;
    const hitCellFormatFn = (col, val) => {
        if (col === "p-value" || col === "Empirical p-value" || col === "FDR"
            || col === "FDR (recruited family)" || col === "padj") return _fmtSci(val, 3);
        if (col === "Propagated score") {
            if (val === "") return "";
            const num = Number(val);
            return isFinite(num) ? num.toExponential(2) : val;
        }
        if (col === "Z-score") return _fmtFixed(val, 2);
        if (col === "log2FC" || col === "Seed score") return _fmtFixed(val, 3);
        return val;
    };
    // A seed clearing FDR is largely CIRCULAR (see the caveat above): it
    // receives restart mass directly and is compared against permutations in
    // which it is usually not a seed. Defaulting this ON when there IS a
    // recruited set puts the 37-ish genes that are the actual discovery in
    // front of the reader, instead of buried under hundreds of expected seed
    // rows - the full table (seeds included) is always one click away, and
    // the CSV download always carries everything regardless of this toggle.
    let hitShowRecruitedOnly = (prop.n_recruited || 0) > 0;
    const hitFilterRow = el("div", { class: "row", style: { alignItems: "center", margin: "0.5rem 0" } });
    const hitFilterLabel = el("label", { style: { display: "flex", alignItems: "center", gap: "0.4rem" } });
    const hitFilterCheckbox = el("input", { type: "checkbox" });
    hitFilterCheckbox.checked = hitShowRecruitedOnly;
    hitFilterLabel.appendChild(hitFilterCheckbox);
    hitFilterLabel.appendChild(document.createTextNode(
        " Show RECRUITED genes only (" + (prop.n_recruited || 0) + " of " + (prop.n_hits || 0) +
        " significant proteins - genes you never uploaded)"));
    hitFilterRow.appendChild(hitFilterLabel);
    hitPane.appendChild(hitFilterRow);
    const hitTableHost = el("div", {});
    hitPane.appendChild(hitTableHost);
    function renderHitTable() {
        hitTableHost.innerHTML = "";
        if (!prop.table) return;
        let tablePayload = prop.table;
        if (hitShowRecruitedOnly) {
            // Filters on the RECRUITED-FAMILY FDR (a seed-free BH family),
            // not the all-protein "FDR" column - a seed clearing the
            // all-protein family threshold is not what this view means by
            // "recruited"; see the caveat banner for why the two families
            // can disagree.
            const rows = (prop.table.rows || []).filter(
                (row) => row["Seed"] === false && row["FDR (recruited family)"] !== null
                    && Number(row["FDR (recruited family)"]) < hitAlpha);
            tablePayload = { columns: prop.table.columns, rows };
        }
        hitTableHost.appendChild(buildTable(tablePayload, {
            maxHeight: 520,
            rowClassFn: (row) => {
                const col = hitShowRecruitedOnly ? "FDR (recruited family)" : "FDR";
                return (row[col] !== null && Number(row[col]) < hitAlpha) ? "row-significant" : "";
            },
            cellFormatFn: hitCellFormatFn,
        }));
    }
    hitFilterCheckbox.addEventListener("change", () => {
        hitShowRecruitedOnly = hitFilterCheckbox.checked;
        renderHitTable();
    });
    renderHitTable();
    const hitDl = el("button", { class: "btn" }, "Download hit table (CSV)");
    hitDl.addEventListener("click", () => downloadFile(
        "/api/transcriptomics/network/download?which=hits", "transcriptomics_network_hits.csv"));
    hitPane.appendChild(hitDl);
    hitPane.appendChild(el("div", { class: "caption" },
        "The download always includes every protein, seeds and recruited alike, regardless of this toggle."));

    // ----- Figures -----
    const figPane = nav.pane("Figures");
    const figState = { nodesShown: 50, fontSize: 12, colorScaleName: ORA_DOT_COLOR_SCALES[0][1], minEdgeConfidence: 0 };
    const figToggles = [];
    const figControls = el("div", { class: "chart-toggle-row" });
    function figOption(label, choices, current, onChange) {
        const w = el("div", { class: "chart-option-control" });
        w.appendChild(el("span", { class: "chart-option-label" }, label));
        const sel = el("select", { class: "chart-option-select" });
        choices.forEach(([lbl, val]) => sel.appendChild(el("option", { value: String(val) }, lbl)));
        sel.value = String(current);
        sel.addEventListener("change", () => { onChange(sel.value); figToggles.forEach((t) => t.redraw()); });
        w.appendChild(sel);
        figControls.appendChild(w);
    }
    figOption("Nodes shown", [["25", 25], ["50", 50], ["100", 100], ["200", 200]], figState.nodesShown,
        (v) => { figState.nodesShown = Number(v); });
    figOption("Font size", [["10", 10], ["12", 12], ["14", 14]], figState.fontSize,
        (v) => { figState.fontSize = Number(v); });
    figOption("Colour scale", ORA_DOT_COLOR_SCALES, figState.colorScaleName,
        (v) => { figState.colorScaleName = v; });
    figOption("Min edge confidence", [["any", 0], ["0.7", 0.7], ["0.9", 0.9]], figState.minEdgeConfidence,
        (v) => { figState.minEdgeConfidence = Number(v); });
    figPane.appendChild(figControls);

    const netToggle = _netpropChartToggle("Show network figure", "propnet",
        (d) => buildNetpropNetwork(d, payload, "score", { ...figState, dbName: DB }),
        _oraChartId("propnet", DB, "score"), DB);
    const histToggle = _netpropChartToggle("Show null distribution", "histogram",
        (d) => buildNetpropHistogram(d, payload, { fontSize: figState.fontSize, dbName: DB }),
        _oraChartId("histogram", DB, "null"), DB);
    figToggles.push(netToggle, histToggle);
    figPane.appendChild(netToggle.el);
    figPane.appendChild(histToggle.el);
    if (nd.skew_caveat) figPane.appendChild(alertBox("warning", nd.skew_caveat));

    // ----- Modules -----
    if (payload.modules && (payload.modules.table.rows || []).length) {
        const mods = payload.modules;
        const modPane = nav.pane("Modules");
        modPane.appendChild(el("div", { class: "caption" },
            "Densely-connected groups among the significant genes. \"Up\"/\"Down\" count how many members " +
            "moved in each direction in YOUR data, so a module mixing both is not a coherently regulated " +
            "unit. Recruited members have no direction, since they were not in your table."));
        if (mods.cap_info && mods.cap_info.capped) {
            modPane.appendChild(alertBox("info",
                "Clustering used the " + mods.cap_info.n_used_for_modules + " strongest of " +
                mods.cap_info.n_hits + " significant genes."));
        }
        modPane.appendChild(buildTable(mods.table, {
            maxHeight: 460,
            cellFormatFn: (col, val) => {
                if (col === "Top enriched FDR" || col === "Top enriched FDR (raw, uncorrected)") return _fmtSci(val, 3);
                if (col === "Mean log2FC") return _fmtFixed(val, 3);
                return val;
            },
        }));
        const modDl = el("button", { class: "btn" }, "Download modules (CSV)");
        modDl.addEventListener("click", () => downloadFile(
            "/api/transcriptomics/network/download?which=modules", "transcriptomics_network_modules.csv"));
        modPane.appendChild(modDl);

        const modNetToggle = _netpropChartToggle("Show module-coloured network", "modulenet",
            (d) => buildNetpropNetwork(d, payload, "module", { ...figState, dbName: DB }),
            _oraChartId("modulenet", DB, "module"), DB);
        modPane.appendChild(modNetToggle.el);

        // ----- Per-module enrichment -----
        const withEnr = (mods.table.rows || [])
            .map((r) => r["Module ID"])
            .filter((mid) => mods.enrichment && mods.enrichment[mid] && (mods.enrichment[mid].table.rows || []).length);
        if (withEnr.length) {
            const enrPane = nav.pane("Module enrichment");
            enrPane.appendChild(el("div", { class: "caption" },
                mods.selection_caveat || ""));
            const sel = el("select", { class: "text-input" });
            withEnr.forEach((mid) => sel.appendChild(el("option", { value: mid }, mid)));
            enrPane.appendChild(el("div", { class: "row" }, el("div", { class: "col" }, sel)));
            const host = el("div", {});
            enrPane.appendChild(host);
            function renderEnr(mid) {
                host.innerHTML = "";
                const e = mods.enrichment[mid];
                if (!e) return;
                host.appendChild(el("div", { class: "caption" }, "Annotation source: " + e.namespace));
                host.appendChild(buildTable(e.table, {
                    maxHeight: 380,
                    cellFormatFn: (col, val) => {
                        if (col === "p-value" || col === "FDR" || col === "Bonferroni") return _fmtSci(val, 3);
                        if (col === "Fold Enrichment") return _fmtFixed(val, 2);
                        return val;
                    },
                }));
                const dl = el("button", { class: "btn" }, "Download this module's enrichment (CSV)");
                dl.addEventListener("click", () => downloadFile(
                    "/api/transcriptomics/network/download?which=enrichment&key=" + encodeURIComponent(mid),
                    "transcriptomics_network_enrichment_" + mid.toLowerCase() + ".csv"));
                host.appendChild(dl);
                const dot = _netpropChartToggle("Show dot plot", "dot",
                    (d) => buildDotPlot(d, e.table, DB, "enrichment", ORA_DOT_COLOR_SCALES[0][1], { topN: 20, fontSize: 12 }),
                    _oraChartId("dot", DB, "enrichment"), DB);
                const bar = _netpropChartToggle("Show bar plot", "bar",
                    (d) => buildBarPlot(d, e.table, DB, "enrichment", ORA_BAR_COLORS[0][1], { topN: 20, fontSize: 12 }),
                    _oraChartId("bar", DB, "enrichment"), DB);
                host.appendChild(dot.el);
                host.appendChild(bar.el);
            }
            sel.addEventListener("change", () => renderEnr(sel.value));
            renderEnr(withEnr[0]);
        }
    }

    // ----- Parameters -----
    const paramPane = nav.pane("Parameters used");
    const paramBox = el("div", { class: "rank-card ora-db-card" });
    [
        ["STRING version", params.string_version],
        ["Confidence cutoff", params.confidence],
        ["Seed selection", "padj <= " + params.padj_cutoff + ", |log2FC| >= " + params.lfc_cutoff + ", " + params.direction],
        ["Seed weighting", params.weighting + (params.ranking_metric ? " (" + params.ranking_metric + ")" : "")],
        ["Restart probability", params.restart_prob],
        ["Permutations", params.n_permutations + " (" + params.permutation_type + ", seed " + params.seed_rng + ")"],
        ["Significance", "FDR < " + params.fdr_alpha],
        ["MCL inflation", params.mcl_inflation],
        ["Min module size", params.min_module_size],
        ["Max hits clustered", params.max_hits_for_modules],
        ["Module labels", params.enrichment_database + ((params.msigdb_collections || []).length ? " (" + params.msigdb_collections.join(", ") + ")" : "")],
    ].forEach(([k, v]) => paramBox.appendChild(el("div", { class: "caption" }, "- " + k + ": " + v)));
    paramPane.appendChild(paramBox);
    paramPane.appendChild(el("div", { class: "caption" },
        "STRING's combined score includes text-mining and cross-organism transfer as well as experimental " +
        "evidence, so modules lean toward well-studied genes. Module labels are descriptive: modules are " +
        "chosen for being connected and carrying signal and only then tested, so their FDR is not a " +
        "hypothesis test."));

    nav.finish();

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.transcriptomics.networkConfig = null;
        STATE.transcriptomics.networkResults = null;
        STATE.transcriptomics.networkRunLocked = false;
        STATE.transcriptomics._networkDraft = null;
        STATE.aaPage = "transcriptomics_network_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToTranscriptomicsActions());
    return wrap;
}


// ---------------------------------------------------------------------------
// GSEA chart builders. Column names differ from ORA's (NES/ES/"FDR q-value"/
// "Leading Edge Size"/"Leading Edge Genes" rather than k/K/n/N/"Fold
// Enrichment"/"FDR"/"Overlapping Genes"), so these are separate builders -
// same pattern as the Burden Test clones - but reuse every generic piece
// (_newOraPlot, _oraChartId, _oraCategoricalChartHeight, _wrapOraLabel,
// _forceLayout, _plotNetwork, _makeNetworkDraggable) so charts look and
// export identically to every other module's.
// ---------------------------------------------------------------------------

function _gseaTopRows(table, n, alpha) {
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : GSEA_DEFAULT_SIG_THRESHOLD;
    const rows = _tableToObjects(table).filter(r =>
        r["FDR q-value"] !== null && r["FDR q-value"] !== "" && isFinite(Number(r["FDR q-value"])) && Number(r["FDR q-value"]) < sigAlpha
    );
    // rows already arrive FDR q-value ascending (then |NES| descending) from the server
    return n === Infinity ? rows : rows.slice(0, n);
}

// GSEA analog of _oraGeneSet/_burdenGeneSet: the leading-edge gene list lives
// in the "Leading Edge Genes" column, comma-joined the same way.
function _gseaGeneSet(row) {
    const raw = String(row["Leading Edge Genes"] || "").trim();
    if (!raw) return new Set();
    return new Set(raw.split(",").map(g => g.trim()).filter(Boolean));
}

const GSEA_DOTBAR_HOVER = "%{customdata[0]}<br>NES: %{customdata[1]:.2f}<br>FDR q: %{customdata[2]:.3e}<br>Leading edge: %{customdata[3]}<extra></extra>";

function buildGseaDotPlot(div, table, slug, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const rows = _gseaTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant gene sets to plot."; return; }
    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const fdrVals = rows.map(r => Number(r["FDR q-value"]));
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const trace = {
        x: rows.map(r => Number(r.NES)),
        y: wrappedNames,
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.NES), Number(r["FDR q-value"]), Number(r["Leading Edge Size"]),
        ]),
        mode: "markers",
        type: "scatter",
        marker: {
            size: _magnitudeMarkerSizes(rows.map(r => Number(r["Leading Edge Size"])), 8, 30),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR q)", thickness: 14, len: 0.7 },
        },
        hovertemplate: GSEA_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "NES", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 105, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "dot", "Dot plot");
}

function buildGseaBarPlot(div, table, slug, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const rows = _gseaTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant gene sets to plot."; return; }
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const nesVals = rows.map(r => Number(r.NES));
    const trace = {
        x: nesVals,
        y: wrappedNames,
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.NES), Number(r["FDR q-value"]), Number(r["Leading Edge Size"]),
        ]),
        type: "bar",
        orientation: "h",
        // Colored by NES sign (up/down), matching how GSEA results are
        // conventionally read, rather than the ORA bar's single flat color.
        marker: { color: nesVals.map(v => v >= 0 ? "#B84A52" : "#1B639C") },
        hovertemplate: GSEA_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "NES", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 70, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "bar", "Bar plot");
}

function buildGseaVolcanoPlot(div, table, slug, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : GSEA_DEFAULT_SIG_THRESHOLD;
    const rows = _tableToObjects(table).filter(r => r.NES !== null && isFinite(Number(r.NES)));
    if (!rows.length) { div.textContent = "No gene sets to plot."; return; }
    const x = rows.map(r => Number(r.NES));
    const y = rows.map(r => -Math.log10(Number(r["FDR q-value"]) + 1e-300));
    const colors = rows.map(r => Number(r["FDR q-value"]) < sigAlpha ? "#c0392b" : "#95a5a6");
    const trace = {
        x, y, mode: "markers", type: "scatter",
        text: rows.map(r => r["Term Name"]),
        marker: { color: colors, size: 8, opacity: 0.8 },
        hovertemplate: "%{text}<br>NES=%{x:.2f}<br>-log10(FDR q)=%{y:.2f}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "NES", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { title: { text: "-log10(FDR q-value)", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: false },
        height: 500,
        margin: { l: 95, r: 110, t: 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "volcano", "Volcano");
}

function buildGseaEnrichmentMap(div, table, slug, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const edgeThreshold = opts.edgeThreshold || 1;
    const rows = _gseaTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant gene sets to plot."; return; }
    if (rows.length < 2) { div.textContent = "Need at least 2 gene sets to draw an enrichment map."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const geneSets = rows.map(_gseaGeneSet);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r["FDR q-value"]));
    const sizeVals = rows.map(r => Number(r["Leading Edge Size"]));
    const nesVals = rows.map(r => Number(r.NES));

    const edges = [];
    for (let i = 0; i < rows.length; i++) {
        for (let j = i + 1; j < rows.length; j++) {
            let overlap = 0;
            geneSets[i].forEach((g) => { if (geneSets[j].has(g)) overlap++; });
            if (overlap < edgeThreshold) continue;
            const union = geneSets[i].size + geneSets[j].size - overlap;
            const jaccard = union > 0 ? overlap / union : 0;
            edges.push({ source: termIds[i], target: termIds[j], weight: jaccard, overlap });
        }
    }

    const pos = _forceLayout(termIds, edges, { iterations: 150 });
    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#B7C6D6", width: 1.5 }, hoverinfo: "skip", showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const nodeTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], nesVals[i], fdrVals[i], sizeVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(sizeVals, 10, 34),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR q)", thickness: 14, len: 0.7 },
            line: { color: "#2c3e50", width: 1 },
        },
        hovertemplate: "%{customdata[0]}<br>NES: %{customdata[1]:.2f}<br>FDR q: %{customdata[2]:.3e}<br>Leading edge: %{customdata[3]}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, nodeTrace], {}, slug, null, "enrichmap", "Enrichment map",
        () => _makeNetworkDraggable(div, [{ traceIndex: 1, nodeIds: termIds }], pos, 0, edges));
}

function buildGseaGeneConceptNetwork(div, table, slug, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    // Gene-concept networks get crowded fast (every gene becomes its own node),
    // so cap top-N tighter than the dot/bar default even if a larger value is set.
    const topN = opts.topN === Infinity ? 15 : Math.min(opts.topN || 20, 15);
    const rows = _gseaTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant gene sets to plot."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === opts.dotColorScale) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r["FDR q-value"]));
    const sizeVals = rows.map(r => Number(r["Leading Edge Size"]));
    const nesVals = rows.map(r => Number(r.NES));
    const geneSets = rows.map(_gseaGeneSet);

    const geneIndex = {};
    const geneIds = [];
    geneSets.forEach((set) => {
        set.forEach((g) => {
            if (!(g in geneIndex)) { geneIndex[g] = "G" + geneIds.length; geneIds.push(g); }
        });
    });

    const nodeIds = termIds.concat(geneIds.map(g => geneIndex[g]));
    const edges = [];
    rows.forEach((r, i) => {
        geneSets[i].forEach((g) => edges.push({ source: termIds[i], target: geneIndex[g], weight: 1 }));
    });

    const pos = _forceLayout(nodeIds, edges, { iterations: 150 });
    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#D8E1EA", width: 1 }, hoverinfo: "skip", showlegend: false,
    };

    const geneTrace = {
        x: geneIds.map(g => pos[geneIndex[g]].x),
        y: geneIds.map(g => pos[geneIndex[g]].y),
        mode: "markers",
        type: "scatter",
        text: geneIds,
        hovertemplate: "%{text}<extra></extra>",
        marker: { size: 8, color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } },
        showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const termTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], nesVals[i], fdrVals[i], sizeVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(sizeVals, 14, 38),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR q)", thickness: 14, len: 0.7 },
            line: { color: "#2c3e50", width: 1.5 },
            symbol: "diamond",
        },
        hovertemplate: "%{customdata[0]}<br>NES: %{customdata[1]:.2f}<br>FDR q: %{customdata[2]:.3e}<br>Leading edge: %{customdata[3]}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, geneTrace, termTrace], {}, slug, null, "network", "Gene-concept network",
        () => _makeNetworkDraggable(
            div,
            [{ traceIndex: 1, nodeIds: geneIds.map(g => geneIndex[g]) }, { traceIndex: 2, nodeIds: termIds }],
            pos, 0, edges
        ));
}

// The signature 3-panel GSEA figure for ONE gene set: running enrichment
// score (top), a rug of where this gene set's members fall in the ranked
// list (middle), and the ranking metric across the whole list (bottom) - all
// sharing one x-axis (gene rank). No other chart in this app uses stacked
// y-axis domains; this is the first one, built directly on Plotly's
// layout.yaxis/yaxis2/yaxis3 domain mechanism (already available - the full
// Plotly bundle is loaded - just never previously used here).
function buildGseaEnrichmentPlot(div, curveData, rankingProfile, termRow, fontSize, slug) {
    if (typeof Plotly === "undefined") return;
    if (!curveData || !rankingProfile) { div.textContent = "No stored curve for this gene set."; return; }
    fontSize = fontSize || 12;
    const ranks = curveData.ranks;
    const values = curveData.values;
    const rankAtMax = curveData.rank_at_max;
    const hitRanks = curveData.hit_ranks || [];
    const sign = curveData.sign;
    const peakIdx = ranks.indexOf(rankAtMax);
    const peakY = peakIdx >= 0 ? values[peakIdx] : 0;
    const termLabel = termRow ? String(termRow["Term Name"] || "") : "";

    const esTrace = {
        x: ranks, y: values, mode: "lines", type: "scatter",
        fill: "tozeroy",
        line: { color: sign >= 0 ? "#B84A52" : "#1B639C", width: 2 },
        fillcolor: sign >= 0 ? "rgba(184,74,82,0.15)" : "rgba(27,99,156,0.15)",
        yaxis: "y",
        hovertemplate: "Rank %{x}<br>Running ES %{y:.3f}<extra></extra>",
        showlegend: false, name: "Running ES",
    };
    const peakTrace = {
        x: [rankAtMax], y: [peakY], mode: "markers", type: "scatter",
        marker: { color: "#222", size: 9, symbol: "diamond" },
        yaxis: "y",
        hovertemplate: "Peak (ES) at rank %{x}<extra></extra>",
        showlegend: false, name: "Peak",
    };
    const rugTrace = {
        x: hitRanks, y: hitRanks.map(() => 0), mode: "markers", type: "scatter",
        marker: { symbol: "line-ns-open", size: 12, color: "#333", line: { width: 1 } },
        yaxis: "y2",
        hovertemplate: "Gene-set member at rank %{x}<extra></extra>",
        showlegend: false, name: "Hits",
    };
    const metricTrace = {
        x: rankingProfile.ranks, y: rankingProfile.values, type: "scatter", mode: "lines",
        fill: "tozeroy", line: { color: "#6C7A89", width: 1 },
        fillcolor: "rgba(108,122,137,0.25)",
        yaxis: "y3",
        hovertemplate: "Rank %{x}<br>Ranking metric %{y:.2f}<extra></extra>",
        showlegend: false, name: "Ranking metric",
    };

    const layout = {
        font: { size: fontSize },
        height: 620,
        margin: { l: 70, r: 30, t: 40, b: 56 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        title: { text: _escapePlotlyText(_wrapOraLabel(termLabel, 90)), font: { size: fontSize + 1 } },
        // anchor:"y3" is required, not cosmetic. Plotly defaults an unanchored
        // xaxis to anchor:"y" - the FIRST y-axis - so the axis line, ticks and
        // title get drawn at the bottom of the ES panel, i.e. straight across
        // the Hits rug beneath it. Anchoring to the bottom panel puts them
        // where a GSEA plot actually wants them: under the whole stack.
        xaxis: { anchor: "y3", title: { text: "Rank in ranked gene list", font: { size: fontSize + 1 } }, tickfont: { size: fontSize } },
        // Three vertically stacked panels sharing the x-axis above, via
        // yaxis/yaxis2/yaxis3 domain bands (all anchor:"x").
        yaxis: { domain: [0.44, 1.0], anchor: "x", title: { text: "Enrichment score", font: { size: fontSize + 1 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis2: { domain: [0.32, 0.40], anchor: "x", range: [-1, 1], fixedrange: true, showticklabels: false, showgrid: false, zeroline: false, title: { text: "Hits", font: { size: fontSize - 1 } } },
        yaxis3: { domain: [0.0, 0.28], anchor: "x", title: { text: "Ranking metric", font: { size: fontSize + 1 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        shapes: [{ type: "line", x0: rankAtMax, x1: rankAtMax, yref: "paper", y0: 0, y1: 1, line: { color: "#888", width: 1, dash: "dot" } }],
    };
    _newOraPlot(div, [esTrace, peakTrace, rugTrace, metricTrace], layout, slug, null, "enrichplot", "Enrichment plot");
}

// ---------------------------------------------------------------------------
// GSEA figures 7-11: ridgeline, multi-set overlay, leading-edge heatmap,
// null-distribution histogram, NES-vs-size diagnostic. Shared helpers first.
// ---------------------------------------------------------------------------

// Metric value for every gene in the ranked list, keyed by symbol. Built once
// per panel render from ranking_full (Transcriptomics_module.run_gsea) -
// NEVER from STATE.transcriptomics.resolved, which knows nothing about a
// background gene list and would silently disagree with the ES curves
// (the server filters the universe BEFORE ranking; the client only ever
// sees a filename/count summary of that background list).
function _gseaMetricBySymbol(rankingFull) {
    const m = new Map();
    if (!rankingFull || !rankingFull.symbols) return m;
    for (let i = 0; i < rankingFull.symbols.length; i++) {
        const s = rankingFull.symbols[i];
        if (!m.has(s)) m.set(s, rankingFull.values[i]);
    }
    return m;
}

// Member metric values for one gene set. "members" mode uses the stored
// curve's hit_ranks (every set member, only available for the top
// GSEA_MAX_CURVES_STORED terms); "edge" mode uses the Leading Edge Genes
// column, which works for every row regardless of whether it has a curve.
function _gseaSetValues(row, curve, rankingFull, metricBySymbol, mode) {
    if (mode === "members" && curve && curve.hit_ranks && rankingFull) {
        return curve.hit_ranks.map(r => rankingFull.values[r]).filter(v => isFinite(v));
    }
    const out = [];
    _gseaGeneSet(row).forEach(g => { const v = metricBySymbol.get(g); if (isFinite(v)) out.push(v); });
    return out;
}

// Gaussian KDE on a fixed grid. Silverman's ROBUST rule (Silverman 1986,
// eq. 3.31): bandwidth = 0.9 * min(sd, IQR/1.34) * n^(-1/5). The 0.9 constant
// is paired specifically with the IQR-robust sigma - 1.06 is the plain
// normal-reference constant (paired with sd alone) and using it here would
// over-smooth by a factor of ~1.06/0.9 (~18%).
function _gseaKde(values, grid) {
    const n = values.length;
    if (!n) return grid.map(() => 0);
    const mean = values.reduce((a, b) => a + b, 0) / n;
    const sd = Math.sqrt(values.reduce((a, b) => a + (b - mean) * (b - mean), 0) / Math.max(n - 1, 1));
    const sorted = values.slice().sort((a, b) => a - b);
    const q = (p) => sorted[Math.min(n - 1, Math.max(0, Math.floor(p * (n - 1))))];
    const iqr = q(0.75) - q(0.25);
    // Explicit branches, not Math.min(sd || Infinity, ...): when every value
    // is identical, sd===0 AND iqr===0, so `sd || Infinity` (0 is falsy in
    // JS) sends BOTH terms to Infinity - and Math.max(sigma, 1e-6) can only
    // enforce a lower bound, it can never pull an Infinity back down. That
    // silently produced h=Infinity and an all-zero density for a degenerate
    // gene set instead of the narrow spike the floor below is meant to give.
    let sigma;
    if (sd > 0 && iqr > 0) sigma = Math.min(sd, iqr / 1.349);
    else if (sd > 0) sigma = sd;
    else if (iqr > 0) sigma = iqr / 1.349;
    else sigma = 1e-6;   // every value identical - no spread at all
    sigma = Math.max(sigma, 1e-6);
    const h = Math.max(0.9 * sigma * Math.pow(n, -0.2), 1e-6);
    const norm = 1 / (n * h * Math.sqrt(2 * Math.PI));
    return grid.map(x => {
        let s = 0;
        for (let i = 0; i < n; i++) { const u = (x - values[i]) / h; s += Math.exp(-0.5 * u * u); }
        return s * norm;
    });
}

const GSEA_DIVERGING_SCALE = [[0, "#1B639C"], [0.5, "#F7F7F7"], [1, "#B84A52"]];
const GSEA_OVERLAY_COLORS = ["#B84A52", "#1B639C", "#4B8B3B", "#8E5FA8", "#C97D1E"];

// ---------------------------------------------------------------------------
// Figure 7 - Ridgeline plot (clusterProfiler ridgeplot equivalent). A hand-
// computed KDE drawn as closed scatter polygons, NOT Plotly's violin trace:
// a violin is anchored symmetrically inside its own category slot and has no
// way to spill into a neighbouring slot, which is the entire defining
// property of a ridgeline (each density overlapping the one below it).
// ---------------------------------------------------------------------------
function buildGseaRidgePlot(div, table, curves, rankingFull, slug, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const mode = opts.population === "members" ? "members" : "edge";
    let rows = _gseaTopRows(table, opts.topN || 20, opts.alpha);
    if (mode === "members") rows = rows.filter(r => (curves || {})[r["Term ID"]]);
    if (!rows.length) {
        div.textContent = mode === "members"
            ? "No significant gene sets with a stored curve to plot (try 'Leading edge only')."
            : "No significant gene sets to plot.";
        return;
    }
    if (rows.length > 30) rows = rows.slice(0, 30);   // readability cap - ridges overlap

    const metricBySymbol = _gseaMetricBySymbol(rankingFull);
    const series = rows.map(r => ({
        row: r,
        vals: _gseaSetValues(r, (curves || {})[r["Term ID"]], rankingFull, metricBySymbol, mode),
    })).filter(s => s.vals.length >= 3);
    if (!series.length) { div.textContent = "Not enough member genes per set to estimate a density."; return; }

    // Ridges ordered by NES: most positive at the TOP (drawn last = on top).
    series.sort((a, b) => Number(a.row.NES) - Number(b.row.NES));

    const all = series.reduce((acc, s) => acc.concat(s.vals), []);
    const lo = Math.min(...all), hi = Math.max(...all);
    const pad = (hi - lo) * 0.08 || 1;
    const G = 200;
    const grid = Array.from({ length: G }, (_, i) => lo - pad + (i * (hi - lo + 2 * pad)) / (G - 1));

    const dens = series.map(s => _gseaKde(s.vals, grid));
    const peak = Math.max(...dens.map(d => Math.max(...d))) || 1;
    const SPACING = 1;
    const AMP = 1.9 * SPACING;   // > SPACING => ridges overlap, which is the point

    const traces = series.map((s, k) => {
        const base = k * SPACING;
        const nes = Number(s.row.NES);
        const col = nes >= 0 ? "#B84A52" : "#1B639C";
        const rgba = nes >= 0 ? "rgba(184,74,82,0.55)" : "rgba(27,99,156,0.55)";
        const xs = grid.concat(grid.slice().reverse());
        const ys = dens[k].map(d => base + (d / peak) * AMP).concat(grid.map(() => base));
        return {
            x: xs, y: ys, type: "scatter", mode: "lines",
            fill: "toself", fillcolor: rgba,
            line: { color: col, width: 1.2, shape: "spline", smoothing: 0.6 },
            name: String(s.row["Term Name"] || ""),
            hovertemplate: _escapePlotlyText(String(s.row["Term Name"] || "")) +
                "<br>NES " + nes.toFixed(2) + "<br>" + s.vals.length + " genes<extra></extra>",
            showlegend: false,
        };
    });

    // Colour key. The ridges themselves carry showlegend:false (one legend
    // entry per gene set would be a 30-row wall restating the y-axis), so the
    // red/blue split is otherwise unexplained on screen - these two zero-point
    // dummy traces are the only thing telling the reader that hue means the
    // sign of NES and nothing else. Only show a side that actually occurs.
    const anyPos = series.some(s => Number(s.row.NES) >= 0);
    const anyNeg = series.some(s => Number(s.row.NES) < 0);
    const keyTrace = (name, col) => ({
        x: [null], y: [null], type: "scatter", mode: "lines",
        line: { color: col, width: 6 }, name: name, showlegend: true, hoverinfo: "skip",
    });
    if (anyPos) traces.push(keyTrace("NES > 0 (enriched toward the top of the ranking)", "#B84A52"));
    if (anyNeg) traces.push(keyTrace("NES < 0 (enriched toward the bottom)", "#1B639C"));

    const labels = series.map(s => _wrapOraLabel(String(s.row["Term Name"] || ""), 44));
    const perRidge = Math.max(30, Math.round(fontSize * 2.6));
    const layout = {
        font: { size: fontSize },
        height: Math.max(420, series.length * perRidge + 190),
        margin: { l: 385, r: 40, t: 26, b: 96 },
        legend: { orientation: "h", x: 0, y: -0.06, yanchor: "top", font: { size: Math.max(9, fontSize - 1) } },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        xaxis: { title: { text: "Ranking metric", font: { size: fontSize + 2 } }, tickfont: { size: fontSize },
                 zeroline: true, zerolinecolor: "#999", zerolinewidth: 1 },
        yaxis: { tickmode: "array", tickvals: series.map((_, k) => k * SPACING), ticktext: labels,
                 tickfont: { size: fontSize }, automargin: true, showgrid: false, zeroline: false,
                 range: [-0.4 * SPACING, (series.length - 1) * SPACING + AMP + 0.3] },
        shapes: [{ type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1,
                   line: { color: "#999", width: 1, dash: "dot" } }],
    };
    _newOraPlot(div, traces, layout, slug, null, "ridgeline", "Ridgeline plot");
}

// ---------------------------------------------------------------------------
// Figure 8 - Multi-set enrichment overlay (gseaplot2-equivalent): 2-5 gene
// sets' running-ES curves on one figure. One shared rug axis with a fixed
// integer range [0, N] and set k's rug at y = N-k-0.5, rather than N axes
// that appear/disappear with the selection - that would make domain
// arithmetic and the export-height contract a function of N.
// ---------------------------------------------------------------------------
function buildGseaMultiEnrichmentPlot(div, curves, rankingProfile, rowsByTermId, termIds, fontSize, slug) {
    if (typeof Plotly === "undefined") return;
    fontSize = fontSize || 12;
    const ids = (termIds || []).filter(t => curves && curves[t]).slice(0, 5);
    if (ids.length < 2) { div.textContent = "Select 2-5 gene sets (with a stored curve) to overlay."; return; }
    const N = ids.length;

    const traces = [];
    ids.forEach((tid, k) => {
        const c = curves[tid];
        const col = GSEA_OVERLAY_COLORS[k % GSEA_OVERLAY_COLORS.length];
        const row = rowsByTermId[tid] || {};
        const label = _wrapOraLabel(String(row["Term Name"] || tid), 40);
        traces.push({
            x: c.ranks, y: c.values, type: "scatter", mode: "lines", yaxis: "y",
            line: { color: col, width: 2 },
            name: label + " (NES " + Number(row.NES).toFixed(2) + ")",
            hovertemplate: "Rank %{x}<br>Running ES %{y:.3f}<extra>" + _escapePlotlyText(String(row["Term Name"] || tid)) + "</extra>",
            showlegend: true,
        });
        const lane = N - k - 0.5;
        traces.push({
            x: c.hit_ranks || [], y: (c.hit_ranks || []).map(() => lane),
            type: "scatter", mode: "markers", yaxis: "y2",
            marker: { symbol: "line-ns-open", size: 11, color: col, line: { width: 1 } },
            hovertemplate: "Member of " + _escapePlotlyText(String(row["Term Name"] || tid)) + " at rank %{x}<extra></extra>",
            showlegend: false,
        });
    });
    traces.push({
        x: rankingProfile.ranks, y: rankingProfile.values, type: "scatter", mode: "lines",
        fill: "tozeroy", line: { color: "#6C7A89", width: 1 }, fillcolor: "rgba(108,122,137,0.25)",
        yaxis: "y3", hovertemplate: "Rank %{x}<br>Ranking metric %{y:.2f}<extra></extra>",
        showlegend: false, name: "Ranking metric",
    });

    const rugH = 0.05 * N;   // 0.10 (N=2) .. 0.25 (N=5)
    const layout = {
        font: { size: fontSize },
        height: 520 + 26 * N,
        margin: { l: 70, r: 30, t: 30, b: 130 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        // anchor:"y3" - see buildGseaEnrichmentPlot: without it the ticks and
        // title render across the Hits rug instead of below the bottom panel.
        xaxis: { anchor: "y3", title: { text: "Rank in ranked gene list", font: { size: fontSize + 1 } }, tickfont: { size: fontSize } },
        yaxis:  { domain: [0.30 + rugH, 1.0], anchor: "x", title: { text: "Running enrichment score", font: { size: fontSize + 1 } },
                  tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis2: { domain: [0.26, 0.26 + rugH], anchor: "x", range: [0, N], fixedrange: true,
                  showticklabels: false, showgrid: false, zeroline: false,
                  title: { text: "Hits", font: { size: fontSize - 1 } } },
        yaxis3: { domain: [0.0, 0.22], anchor: "x", title: { text: "Ranking metric", font: { size: fontSize + 1 } },
                  tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        legend: { orientation: "h", x: 0, y: -0.16, font: { size: Math.max(9, fontSize - 1) } },
        shapes: ids.map((tid, k) => ({
            type: "line", x0: curves[tid].rank_at_max, x1: curves[tid].rank_at_max,
            yref: "y domain", y0: 0, y1: 1,
            line: { color: GSEA_OVERLAY_COLORS[k % GSEA_OVERLAY_COLORS.length], width: 1, dash: "dot" },
        })),
    };
    _newOraPlot(div, traces, layout, slug, null, "multienrich", "Multi-set enrichment plot");
}

// ---------------------------------------------------------------------------
// Figure 9 - Leading-edge heatmap (enrichplot heatplot equivalent). First
// `heatmap` trace in the app. Rows capped at 25 sets; columns capped
// (default 60) chosen by membership-count descending (keeps the genes that
// actually connect sets, the interesting part of a heatplot) or by strongest
// metric magnitude. No hierarchical clustering: that needs a distance
// matrix + linkage (untested JS, or an unsafe scipy call) and nowhere to put
// a dendrogram in this layout - deterministic sort order instead.
// ---------------------------------------------------------------------------
function buildGseaLeadingEdgeHeatmap(div, table, rankingFull, slug, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const maxGenes = opts.maxGenes || 60;
    let rows = _gseaTopRows(table, Math.min(opts.topN || 20, 25), opts.alpha);
    if (!rows.length) { div.textContent = "No significant gene sets to plot."; return; }

    const metricBySymbol = _gseaMetricBySymbol(rankingFull);
    if (!metricBySymbol.size) { div.textContent = "No ranking data available for this run."; return; }
    const sets = rows.map(r => ({ row: r, genes: _gseaGeneSet(r) }));

    const count = new Map();
    sets.forEach(s => s.genes.forEach(g => { if (metricBySymbol.has(g)) count.set(g, (count.get(g) || 0) + 1); }));
    if (!count.size) { div.textContent = "No leading-edge genes could be matched to the ranked list."; return; }

    let genes = Array.from(count.keys());
    if (opts.geneMode === "metric") {
        genes.sort((a, b) => Math.abs(metricBySymbol.get(b)) - Math.abs(metricBySymbol.get(a)));
    } else {
        genes.sort((a, b) => (count.get(b) - count.get(a)) ||
                             (Math.abs(metricBySymbol.get(b)) - Math.abs(metricBySymbol.get(a))));
    }
    genes = genes.slice(0, maxGenes);
    genes.sort((a, b) => metricBySymbol.get(b) - metricBySymbol.get(a));   // columns: metric descending

    const ordered = opts.rowOrder === "fdr" ? sets : sets.slice().sort((a, b) => Number(b.row.NES) - Number(a.row.NES));
    const yLabels = ordered.map(s => _wrapOraLabel(String(s.row["Term Name"] || ""), 44));
    const z = ordered.map(s => genes.map(g => (s.genes.has(g) ? metricBySymbol.get(g) : null)));
    const lim = Math.max(...genes.map(g => Math.abs(metricBySymbol.get(g)))) || 1;

    const trace = {
        type: "heatmap", x: genes, y: yLabels, z: z,
        colorscale: GSEA_DIVERGING_SCALE, zmid: 0, zmin: -lim, zmax: lim,
        xgap: 1, ygap: 1, hoverongaps: false,
        colorbar: { title: "Ranking<br>metric", thickness: 14, len: 0.7 },
        hovertemplate: "%{y}<br>%{x}<br>Ranking metric: %{z:.3f}<extra></extra>",
    };
    const perRow = Math.max(24, Math.round(fontSize * 1.9));
    const layout = {
        font: { size: fontSize },
        height: Math.max(420, ordered.length * perRow + 210),
        margin: { l: 385, r: 110, t: 24, b: 160 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        xaxis: { tickangle: -60, tickfont: { size: Math.max(8, fontSize - 2) }, automargin: true,
                 side: "bottom", title: { text: "Leading-edge genes (ranked-metric order)", font: { size: fontSize + 1 } } },
        yaxis: { tickfont: { size: fontSize }, automargin: true, autorange: "reversed" },
    };
    _newOraPlot(div, [trace], layout, slug, null, "heatplot", "Leading-edge heatmap");
}

// ---------------------------------------------------------------------------
// Figure 10 - Null-distribution histogram. Mirrors buildNetpropHistogram
// (histogram + shapes line + annotations), sharing the selected gene set
// with the enrichment plot.
// ---------------------------------------------------------------------------
function buildGseaNullHistogram(div, block, termRow, fontSize, slug) {
    if (typeof Plotly === "undefined") return;
    fontSize = fontSize || 12;
    if (!termRow) { div.textContent = "Select a gene set."; return; }
    // JSON object keys are strings - the server writes str(size) deliberately.
    const samples = ((block && block.null_nes_by_size) || {})[String(termRow["Size"])];
    if (!samples || !samples.length) { div.textContent = "No stored permutation null for this gene set."; return; }

    const nes = Number(termRow.NES);
    const col = nes >= 0 ? "#B84A52" : "#1B639C";
    const trace = {
        x: samples, type: "histogram", nbinsx: 60,
        marker: { color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } },
        hovertemplate: "Null NES %{x:.2f}<br>%{y} permutations<extra></extra>",
    };
    // Plotly's own autorange DOES stretch to cover a shape line past the
    // histogram's data range, but with near-zero padding when that shape
    // sits exactly at the new extreme - pinning the observed-NES line onto
    // the axis border itself, indistinguishable from it. Confirmed live:
    // an NES more extreme than the whole null (exactly the significant
    // case this line exists to highlight) rendered invisibly this way.
    // An explicit range with real padding keeps it visibly clear of the
    // edge regardless of how extreme the observed NES is.
    const sampleMin = Math.min(...samples);
    const sampleMax = Math.max(...samples);
    const rangeMin = Math.min(sampleMin, nes, 0);
    const rangeMax = Math.max(sampleMax, nes, 0);
    const rangePad = (rangeMax - rangeMin || 1) * 0.06;
    const layout = {
        font: { size: fontSize },
        height: 460, margin: { l: 70, r: 30, t: 40, b: 56 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        title: { text: _escapePlotlyText(_wrapOraLabel(String(termRow["Term Name"] || ""), 90)), font: { size: fontSize + 1 } },
        xaxis: { title: { text: "Permutation null NES (gene sets of size " + termRow["Size"] + ")", font: { size: fontSize + 2 } },
                 tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc",
                 range: [rangeMin - rangePad, rangeMax + rangePad] },
        // Linear y, deliberately unlike buildNetpropHistogram's log y: a
        // 500-sample thinned null in 60 bins has counts ~0-15, and a log
        // axis would visually inflate exactly the tails being judged here.
        yaxis: { title: { text: "Permutations", font: { size: fontSize + 2 } }, tickfont: { size: fontSize } },
        shapes: [
            { type: "line", x0: 0, x1: 0, yref: "paper", y0: 0, y1: 1, line: { color: "#999", width: 1, dash: "dot" } },
            { type: "line", x0: nes, x1: nes, yref: "paper", y0: 0, y1: 1, line: { color: col, width: 2 } },
        ],
        annotations: [{
            x: nes, yref: "paper", y: 1, showarrow: false,
            text: "Observed NES " + nes.toFixed(2) + " (p = " + Number(termRow["p-value"]).toExponential(2) + ")",
            font: { size: Math.max(9, fontSize - 2), color: col },
            xanchor: nes >= 0 ? "right" : "left",
        }],
    };
    _newOraPlot(div, [trace], layout, slug, null, "nullhist", "Null distribution");
}

// Equal-count size bins for the NES-vs-size diagnostic below. Requires at
// least minPerBin points per bin (default 8) - below that a "10-90th
// percentile" band is really just min-max on 2-3 points, which mislabels
// noise as a summary. nBins < 2 signals "too few gene sets to bin at all";
// the caller must skip drawing the band/median in that case rather than
// drawing a misleading one.
function _gseaSizeBins(rows, minPerBin) {
    minPerBin = minPerBin || 8;
    const sorted = rows.slice().sort((a, b) => Number(a.Size) - Number(b.Size));
    const n = sorted.length;
    const nBins = Math.min(8, Math.floor(n / minPerBin));
    if (nBins < 2) return { bx: [], p10: [], p50: [], p90: [], nBins: nBins, perBin: [] };
    const pct = (arr, p) => { const s = arr.slice().sort((a, b) => a - b);
        return s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))]; };
    const bx = [], p10 = [], p50 = [], p90 = [], perBin = [];
    for (let b = 0; b < nBins; b++) {
        const chunk = sorted.slice(Math.floor(b * n / nBins), Math.floor((b + 1) * n / nBins));
        if (!chunk.length) continue;
        const v = chunk.map(r => Number(r.NES));
        bx.push(chunk.reduce((a, r) => a + Number(r.Size), 0) / chunk.length);
        p10.push(pct(v, 0.10)); p50.push(pct(v, 0.50)); p90.push(pct(v, 0.90));
        perBin.push(chunk.length);
    }
    return { bx, p10, p50, p90, nBins, perBin };
}

// ---------------------------------------------------------------------------
// Figure 11 - NES vs gene-set-size diagnostic. Zero backend change (Size and
// NES are already in the table). A binned 10-90 percentile band + median
// line is what turns a raw scatter into a readable diagnostic of whether
// normalization actually removed the size dependence - but only when there
// are enough gene sets to bin meaningfully (see _gseaSizeBins).
// ---------------------------------------------------------------------------
function buildGseaNesVsSizePlot(div, table, slug, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const rows = _tableToObjects(table).filter(r => isFinite(Number(r.NES)) && Number(r.Size) > 0);
    if (!rows.length) { div.textContent = "No gene sets to plot."; return; }
    const scale = ORA_DOT_COLOR_SCALES.find(s => s[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];

    const bins = _gseaSizeBins(rows);
    const traces = [];
    if (bins.nBins >= 2) {
        traces.push(
            { x: bins.bx, y: bins.p10, type: "scatter", mode: "lines", line: { width: 0 }, showlegend: false, hoverinfo: "skip" },
            { x: bins.bx, y: bins.p90, type: "scatter", mode: "lines", line: { width: 0 }, fill: "tonexty",
              fillcolor: "rgba(120,140,160,0.18)", name: "10-90th percentile", hoverinfo: "skip" },
            { x: bins.bx, y: bins.p50, type: "scatter", mode: "lines+markers", line: { color: "#33475B", width: 2 },
              marker: { size: 7, color: "#33475B" }, name: "Median NES",
              hovertemplate: "Mean size %{x:.0f}<br>Median NES %{y:.2f}<extra></extra>" },
        );
    }
    traces.push(
        { x: rows.map(r => Number(r.Size)), y: rows.map(r => Number(r.NES)),
          type: "scatter", mode: "markers", name: "Gene sets",
          text: rows.map(r => r["Term Name"]),
          marker: { size: 8, opacity: 0.75,
                    color: rows.map(r => -Math.log10(Number(r["FDR q-value"]) + 1e-300)),
                    colorscale: scale[2], reversescale: scale[3], showscale: true,
                    colorbar: { title: "-log10(FDR q)", thickness: 14, len: 0.7 } },
          hovertemplate: "%{text}<br>Size %{x}<br>NES %{y:.2f}<extra></extra>" },
    );
    const layout = {
        font: { size: fontSize }, height: 520, margin: { l: 70, r: 110, t: 24, b: 60 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
        xaxis: { type: "log", title: { text: "Gene-set size (genes in the ranked list)", font: { size: fontSize + 2 } }, tickfont: { size: fontSize } },
        yaxis: { title: { text: "NES", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#999" },
        legend: { orientation: "h", x: 0, y: -0.14, font: { size: Math.max(9, fontSize - 1) } },
    };
    if (bins.nBins < 2) {
        layout.annotations = [{
            xref: "paper", yref: "paper", x: 0.02, y: 0.98, showarrow: false,
            align: "left", font: { size: Math.max(10, fontSize - 1), color: "#7A5A2E" },
            text: "Too few gene sets (" + rows.length + ") to summarize by size - showing individual points only.",
        }];
    }
    _newOraPlot(div, traces, layout, slug, null, "nessize", "NES vs set size");
}

function renderGenomicsUploadScreen(moduleName) {
    const wrap = el("div", {});
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle" }, "Genomics - " + moduleName));
    wrap.appendChild(el("div", { class: "step-subheader" }, "Upload gene table"));

    const row = el("div", { class: "row" });
    const symBtn = el("button", { class: "btn" }, "Gene Symbol");
    symBtn.addEventListener("click", () => { STATE.genomicsIdentifierType = "gene_symbol"; render(); });
    const ensBtn = el("button", { class: "btn" }, "Ensembl Gene ID");
    ensBtn.addEventListener("click", () => { STATE.genomicsIdentifierType = "ensembl_gene_id"; render(); });
    row.appendChild(el("div", { class: "col" }, symBtn));
    row.appendChild(el("div", { class: "col" }, ensBtn));
    row.appendChild(el("div", { class: "col", style: { flex: "2" } }));
    wrap.appendChild(row);

    if (STATE.genomicsIdentifierType === "gene_symbol") {
        wrap.appendChild(fileInputWidget({ accept: ".csv", buttonLabel: "Choose CSV file" }));
        wrap.appendChild(el("div", { class: "caption" }, "Required column: Gene Symbol. Optional column: UniProt ID"));
    } else if (STATE.genomicsIdentifierType === "ensembl_gene_id") {
        wrap.appendChild(fileInputWidget({ accept: ".csv", buttonLabel: "Choose CSV file" }));
        wrap.appendChild(el("div", { class: "caption" }, "Required column: Ensembl Gene ID. Optional column: UniProt ID"));
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderOraOrganismSelector() {
    const wrap = el("div", {});
    const input = el("input", { class: "text-input", type: "text", placeholder: "Species / organism (name or UniProt taxon ID)" });
    const resultsSelect = el("select", { class: "hidden" });
    const selectedCaption = el("div", { class: "caption" });

    function refreshCaption() {
        selectedCaption.textContent = STATE.ora.taxonId ? "Selected: " + STATE.ora.taxonName + " (taxon " + STATE.ora.taxonId + ")" : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = await apiJson("/api/contexts/search?database=" + encodeURIComponent("Rhea + UniProt") + "&query=" + encodeURIComponent(query) + "&limit=10");
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                (data.results || []).forEach((rec, idx) => {
                    const label = rec.context_id && rec.context_name ? rec.context_id + " | " + rec.context_name : (rec.context_name || rec.context_id || "");
                    resultsSelect.appendChild(el("option", { value: String(idx) }, label));
                });
                resultsSelect._records = data.results || [];
                resultsSelect.classList.toggle("hidden", (data.results || []).length === 0);
            } catch (err) { /* keep prior dropdown contents on transient search errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        STATE.ora.taxonId = String(rec.context_id || "");
        STATE.ora.taxonName = String(rec.context_name || "");
        render();
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function renderOraUploadSection() {
    const wrap = el("div", {});
    const idType = STATE.ora.identifierType;
    let caption, label;
    if (idType === "gene_symbol") { caption = "Required column: Gene Symbol. Optional column: UniProt ID."; label = "Upload gene table"; }
    else if (idType === "ensembl_gene_id") { caption = "Required column: Ensembl Gene ID. Optional column: UniProt ID."; label = "Upload gene table"; }
    else { caption = "Required columns: Sequence, Sequence Type (dna or protein). Optional column: UniProt ID."; label = "Upload sequence table"; }

    wrap.appendChild(el("div", { class: "section-label" }, label));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.ora.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            try {
                let data;
                if (idType === "sequence") {
                    data = await apiUpload("/api/ora/seq/upload", formData);
                    STATE.ora.seqInfo = { nToBlast: data.n_to_blast };
                } else {
                    formData.append("identifier_type", idType);
                    data = await apiUpload("/api/ora/genes/upload", formData);
                }
                STATE.ora.geneTable = data.table;
                STATE.ora.uploadedFilename = file.name;
                STATE.ora.duplicatesRemoved = data.duplicates_removed;
                STATE.ora.resolved = null;
                STATE.ora.finalTable = null;
                STATE.ora.ambiguousChoices = {};
                STATE.ora.blast = { selection: {} };
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));
    wrap.appendChild(el("div", { class: "caption" }, caption));

    if (STATE.ora.geneTable) {
        if (STATE.ora.duplicatesRemoved > 0) {
            const noun = idType === "sequence" ? "sequence(s)" : "gene(s)";
            wrap.appendChild(alertBox("info", "Removed " + STATE.ora.duplicatesRemoved + " duplicate " + noun + "; kept the most informative occurrence of each."));
        }
        const maxHeight = Math.min(400, Math.max(200, (STATE.ora.geneTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.ora.geneTable, {
            maxHeight,
            cellFormatFn: (col, val) => (col === "Sequence" ? truncate(val, 40) : val),
        }));
    }

    return wrap;
}

function renderGeneResolveControls() {
    const wrap = el("div", {});
    const disabled = !STATE.ora.geneTable || (STATE.ora.identifierType === "gene_symbol" && !STATE.ora.taxonId);
    const btn = el("button", { class: "btn" }, "Resolve genes");
    btn.disabled = disabled;
    btn.addEventListener("click", async () => {
        btn.disabled = true;
        const spinner = el("div", { class: "caption" }, "Resolving gene identifiers via UniProt...");
        wrap.appendChild(spinner);
        try {
            const { job_id } = await apiPostJson("/api/ora/genes/resolve", { taxon_id: STATE.ora.taxonId });
            pollJob(job_id, {
                onDone: (result) => {
                    STATE.ora.resolved = result;
                    STATE.ora.ambiguousChoices = {};
                    if (result.auto_final) STATE.ora.finalTable = result.auto_final;
                    render();
                },
                onError: (err) => { alert(err); render(); },
            });
        } catch (err) {
            alert(err.message);
            render();
        }
    });
    wrap.appendChild(btn);
    return wrap;
}

function renderGeneResults() {
    const wrap = el("div", {});
    const resolved = STATE.ora.resolved;
    wrap.appendChild(el("div", { class: "section-label" }, "Resolved genes"));
    wrap.appendChild(buildTable(resolved.resolved_table, {}));
    wrap.appendChild(metricCards([
        ["Resolved", resolved.metrics.resolved],
        ["Resolved (TrEMBL)", resolved.metrics.resolved_unreviewed],
        ["Resolved via NCBI", resolved.metrics.resolved_via_ncbi],
        ["User provided", resolved.metrics.user_provided],
        ["Not found", resolved.metrics.not_found],
    ], 5));

    if (resolved.ambiguous && resolved.ambiguous.length > 0) {
        wrap.appendChild(el("div", { class: "section-label" }, "Ambiguous genes"));
        wrap.appendChild(alertBox("warning", resolved.ambiguous.length + " gene(s) returned multiple UniProt candidates - choose one for each or discard."));

        const discardAllBtn = el("button", { class: "btn" }, "Discard all unresolved");
        discardAllBtn.addEventListener("click", async () => {
            try {
                const data = await apiPostJson("/api/ora/genes/discard-all", {});
                STATE.ora.finalTable = data.final_table;
                render();
            } catch (err) { alert(err.message); }
        });
        wrap.appendChild(discardAllBtn);

        wrap.appendChild(el("div", { style: { fontWeight: "700", margin: "0.8rem 0 0.4rem 0" } }, "Select one UniProt entry per gene, or choose to discard:"));

        resolved.ambiguous.forEach((amb) => {
            const row = el("div", { class: "row", style: { alignItems: "center", margin: "0.3rem 0" } });
            row.appendChild(el("div", { class: "col", style: { flex: "0 0 200px", fontWeight: "700" } }, amb.gene_id));
            const select = el("select", { class: "col" });
            amb.candidates.forEach((c, idx) => {
                select.appendChild(el("option", { value: String(idx) }, c.accession + " - " + c.entry_name + " - " + c.protein_name));
            });
            select.appendChild(el("option", { value: "discard" }, "Discard this value"));
            select.value = STATE.ora.ambiguousChoices[amb.gene_id] !== undefined ? STATE.ora.ambiguousChoices[amb.gene_id] : "discard";
            select.addEventListener("change", () => { STATE.ora.ambiguousChoices[amb.gene_id] = select.value; });
            row.appendChild(select);
            wrap.appendChild(row);
        });

        const applyBtn = el("button", { class: "btn" }, "Apply selections");
        applyBtn.addEventListener("click", async () => {
            const choices = resolved.ambiguous.map((amb) => {
                const choice = STATE.ora.ambiguousChoices[amb.gene_id];
                if (choice === undefined || choice === "discard") return { gene_id: amb.gene_id, discard: true };
                const c = amb.candidates[parseInt(choice, 10)];
                return { gene_id: amb.gene_id, accession: c.accession, entry_name: c.entry_name, protein_name: c.protein_name };
            });
            try {
                const data = await apiPostJson("/api/ora/genes/apply-selections", { choices });
                STATE.ora.finalTable = data.final_table;
                render();
            } catch (err) { alert(err.message); }
        });
        wrap.appendChild(applyBtn);
    }

    if (STATE.ora.finalTable) {
        // The "Final" table only differs from "Resolved genes" when the user had to
        // resolve ambiguity; with no ambiguity it's identical, so don't repeat it.
        const hadAmbiguity = resolved.ambiguous && resolved.ambiguous.length > 0;
        if (hadAmbiguity) {
            wrap.appendChild(el("div", { class: "section-label" }, "Final gene → UniProt table"));
            wrap.appendChild(buildTable(STATE.ora.finalTable, {}));
        }
        const dlBtn = el("button", { class: "btn" }, "Download resolved genes table");
        dlBtn.addEventListener("click", () => downloadFile("/api/ora/download/genes.csv", "genes_resolved.csv"));
        wrap.appendChild(dlBtn);
        wrap.appendChild(renderContinueToAnalysisButton());
    }

    return wrap;
}

function renderContinueToAnalysisButton() {
    // Identifier-type-agnostic handoff point: this only requires finalTable
    // (a confirmed UniProt ID list), regardless of whether it came from Gene
    // Symbol, Ensembl, or BLAST resolution. Routes to whichever analysis
    // (ORA or Gene-set Burden Test) the active genomics module is.
    if (!STATE.ora.finalTable || !STATE.ora.finalTable.rows.length) return null;
    const module = STATE.genomicsModuleName;
    const isBurden = module === "Gene-set Burden Test";
    const isNetprop = module === "Network Propagation / Module Detection";
    const isAnnotation = module === "Gene Function Annotation";
    // Burden requires a confirmed score column; netprop's and annotation's
    // scores/params are optional/have defaults, so neither needs such a gate.
    if (isBurden && !STATE.ora.scoreSummary) return null;
    const pageByModule = isBurden ? "genomics_burden_config"
        : isNetprop ? "genomics_netprop_config"
        : isAnnotation ? "genomics_annotation_config"
        : "genomics_ora_config";
    const wrap = el("div", { style: { marginTop: "0.8rem" } });
    const btn = el("button", { class: "btn btn-primary" }, "Continue to analysis");
    btn.addEventListener("click", () => {
        STATE.aaPage = pageByModule;
        render();
    });
    wrap.appendChild(btn);
    return wrap;
}

function renderBlastControls() {
    const wrap = el("div", {});
    const hasPending = STATE.ora.resolved && STATE.ora.resolved.metrics && STATE.ora.resolved.metrics.blast_pending > 0;

    if (STATE.ora.geneTable && !hasPending) {
        const nToBlast = STATE.ora.seqInfo ? STATE.ora.seqInfo.nToBlast : STATE.ora.geneTable.rows.length;
        wrap.appendChild(alertBox("info", nToBlast + " sequence(s) will be submitted to NCBI BLAST (swissprot database). All jobs are submitted at once, then polled in parallel for up to 3 minutes. If NCBI is slow, jobs are saved and a Resume button lets you collect results without resubmitting."));
    }

    const progressWrap = el("div", { class: "progress-wrap hidden" });
    const progressInner = el("div", { class: "progress-bar-inner", style: { width: "0%" } });
    progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
    const progressLabel = el("div", { class: "progress-label" }, "");
    progressWrap.appendChild(progressLabel);
    wrap.appendChild(progressWrap);

    function runBlast() {
        progressWrap.classList.remove("hidden");
        runBtn.disabled = true;
        if (resumeBtn) resumeBtn.disabled = true;
        apiPostJson("/api/ora/seq/blast", {}).then(({ job_id }) => {
            pollJob(job_id, {
                onProgress: (job) => {
                    progressInner.style.width = Math.round((job.progress || 0) * 100) + "%";
                    progressLabel.textContent = job.message || "";
                },
                onDone: (result) => {
                    STATE.ora.resolved = { table: result.table, metrics: result.metrics };
                    STATE.ora.finalTable = null;
                    STATE.ora.blast = { selection: {} };
                    render();
                },
                onError: (err) => { alert("BLAST failed: " + err); render(); },
            });
        }).catch((err) => { alert(err.message); render(); });
    }

    const runBtn = el("button", { class: "btn" }, "Run homology search (BLAST)");
    runBtn.disabled = !STATE.ora.geneTable || hasPending;
    runBtn.addEventListener("click", runBlast);
    wrap.appendChild(runBtn);

    let resumeBtn = null;
    if (hasPending && STATE.ora.geneTable) {
        resumeBtn = el("button", { class: "btn" }, "Resume - collect finished BLAST jobs");
        resumeBtn.addEventListener("click", runBlast);
        wrap.appendChild(resumeBtn);
    }

    return wrap;
}

function renderBlastResults() {
    const wrap = el("div", {});
    const resolved = STATE.ora.resolved;
    const m = resolved.metrics;
    wrap.appendChild(el("div", { class: "section-label" }, "Homology search results"));
    wrap.appendChild(metricCards([
        ["Resolved (homology)", m.resolved_homology], ["User provided", m.user_provided], ["Not found", m.not_found],
        ["BLAST failed", m.blast_failed], ["Still running", m.blast_pending],
    ], 5));

    if (m.blast_pending > 0) {
        wrap.appendChild(alertBox("info", m.blast_pending + " job(s) are still running on NCBI. Results are saved - click Resume below to collect them when ready."));
    }

    wrap.appendChild(buildTable(resolved.table, {
        rowClassFn: (row) => (row["Confidence"] === "low" ? "row-low-confidence" : null),
        cellFormatFn: (col, val) => (col === "Sequence" ? truncate(val, 40) : val),
    }));

    const dlAllBtn = el("button", { class: "btn" }, "Download all results");
    dlAllBtn.addEventListener("click", () => downloadFile("/api/ora/download/blast_all.csv", "sequences_blast_resolved.csv"));
    wrap.appendChild(dlAllBtn);

    wrap.appendChild(el("hr", { class: "divider" }));
    wrap.appendChild(el("div", { style: { fontWeight: "700" } }, "Select sequences to include in the analysis"));
    wrap.appendChild(el("div", { class: "caption" }, "Resolved rows are pre-checked. Uncheck any you want to exclude. Low-confidence rows (red) are pre-checked but flagged for review."));

    const includeStatuses = { resolved_homology: true, user_provided: true };
    const checksWrap = el("div", {});

    function countSelected() {
        return Object.keys(STATE.ora.blast.selection).filter((k) => STATE.ora.blast.selection[k]).length;
    }

    resolved.table.rows.forEach((row, idx) => {
        const status = row["Resolution_status"] || "";
        const seq = String(row["Sequence"] || "");
        const seqType = row["Sequence Type"] || "";
        const uniprot = row["UniProt ID"] || "";
        const pct = row["Percent Identity"] || "";
        const evalue = row["E-value"] || "";
        const seqDisplay = truncate(seq, 40);
        let labelText;
        if (uniprot) {
            labelText = seqDisplay + " -> " + uniprot + " (" + seqType;
            if (pct) labelText += ", " + pct + "% identity";
            if (evalue) labelText += ", e=" + evalue;
            labelText += ")";
        } else {
            labelText = seqDisplay + " (" + seqType + ") - " + status;
            if (status === "blast_pending") labelText += " (pending)";
        }
        const defaultChecked = !!includeStatuses[status];
        if (STATE.ora.blast.selection[idx] === undefined) STATE.ora.blast.selection[idx] = defaultChecked;

        const rowEl = el("label", { class: "checkbox-row" + (status === "blast_pending" ? " disabled" : "") });
        const cb = el("input", { type: "checkbox" });
        cb.checked = STATE.ora.blast.selection[idx];
        cb.disabled = status === "blast_pending";
        cb.addEventListener("change", () => {
            STATE.ora.blast.selection[idx] = cb.checked;
            confirmBtn.textContent = "Confirm selection (" + countSelected() + " row(s))";
            confirmBtn.disabled = countSelected() === 0;
        });
        rowEl.appendChild(cb);
        rowEl.appendChild(document.createTextNode(" " + labelText));
        checksWrap.appendChild(rowEl);
    });
    wrap.appendChild(checksWrap);

    const confirmBtn = el("button", { class: "btn" }, "Confirm selection (" + countSelected() + " row(s))");
    confirmBtn.disabled = countSelected() === 0;
    confirmBtn.addEventListener("click", async () => {
        const selectedIndices = Object.keys(STATE.ora.blast.selection)
            .filter((k) => STATE.ora.blast.selection[k])
            .map((k) => parseInt(k, 10));
        try {
            const data = await apiPostJson("/api/ora/seq/confirm", { selected_indices: selectedIndices });
            STATE.ora.finalTable = data.final_table;
            render();
        } catch (err) { alert(err.message); }
    });
    wrap.appendChild(confirmBtn);

    if (STATE.ora.finalTable && STATE.ora.finalTable.rows.length) {
        wrap.appendChild(el("hr", { class: "divider" }));
        wrap.appendChild(el("div", { style: { fontWeight: "700" } }, "Confirmed selection for ORA"));
        wrap.appendChild(buildTable(STATE.ora.finalTable, {
            cellFormatFn: (col, val) => (col === "Sequence" ? truncate(val, 40) : val),
        }));
        const dlConfirmedBtn = el("button", { class: "btn" }, "Download confirmed selection");
        dlConfirmedBtn.addEventListener("click", () => downloadFile("/api/ora/download/blast_confirmed.csv", "sequences_blast_confirmed.csv"));
        wrap.appendChild(dlConfirmedBtn);
        wrap.appendChild(renderContinueToAnalysisButton());
    }

    return wrap;
}

function renderGenomicsOraScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Over-Representation Analysis"));

    const typeRow = el("div", { class: "row" });
    function typeBtn(label, type) {
        const btn = el("button", { class: "btn" }, label);
        if (STATE.ora.identifierType === type) btn.style.fontWeight = "800";
        btn.addEventListener("click", () => {
            if (STATE.ora.identifierType !== type) {
                STATE.ora = freshOraState();
                STATE.ora.identifierType = type;
                render();
            }
        });
        return btn;
    }
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Gene Symbol", "gene_symbol")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Ensembl Gene ID", "ensembl_gene_id")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Homology (BLAST)", "sequence")));
    wrap.appendChild(typeRow);

    if (!STATE.ora.identifierType) {
        wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
        wrap.appendChild(backToGenomicsActions());
        return wrap;
    }

    if (STATE.ora.identifierType === "gene_symbol") {
        wrap.appendChild(el("div", { class: "section-label" }, "Select organism"));
        wrap.appendChild(renderOraOrganismSelector());
    }

    wrap.appendChild(renderOraUploadSection());

    if (STATE.ora.identifierType === "sequence") wrap.appendChild(renderBlastControls());
    else wrap.appendChild(renderGeneResolveControls());

    if (STATE.ora.resolved) {
        if (STATE.ora.identifierType === "sequence") wrap.appendChild(renderBlastResults());
        else wrap.appendChild(renderGeneResults());
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderGenomicsBurdenScreen() {
    // Retrieval is fully shared with ORA (same STATE.ora, same organism
    // selector, same resolve/ambiguity flow) - only the upload section adds a
    // required variant-score column, and there is no Homology/BLAST path
    // (sequences aren't a stable identifier to key a per-gene score to).
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene-set Burden Test"));

    const typeRow = el("div", { class: "row" });
    function typeBtn(label, type) {
        const btn = el("button", { class: "btn" }, label);
        if (STATE.ora.identifierType === type) btn.style.fontWeight = "800";
        btn.addEventListener("click", () => {
            if (STATE.ora.identifierType !== type) {
                STATE.ora = freshOraState();
                STATE.ora.identifierType = type;
                render();
            }
        });
        return btn;
    }
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Gene Symbol", "gene_symbol")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Ensembl Gene ID", "ensembl_gene_id")));
    typeRow.appendChild(el("div", { class: "col", style: { flex: "1" } }));
    wrap.appendChild(typeRow);

    if (!STATE.ora.identifierType) {
        wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
        wrap.appendChild(backToGenomicsActions());
        return wrap;
    }

    if (STATE.ora.identifierType === "gene_symbol") {
        wrap.appendChild(el("div", { class: "section-label" }, "Select organism"));
        wrap.appendChild(renderOraOrganismSelector());
    }

    wrap.appendChild(renderBurdenUploadSection());
    wrap.appendChild(renderGeneResolveControls());

    if (STATE.ora.resolved) {
        wrap.appendChild(renderGeneResults());
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderBurdenUploadSection() {
    const wrap = el("div", {});
    const idType = STATE.ora.identifierType;
    const idColLabel = idType === "gene_symbol" ? "Gene Symbol" : "Ensembl Gene ID";

    wrap.appendChild(el("div", { class: "section-label" }, "Upload gene table"));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.ora.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            formData.append("identifier_type", idType);
            try {
                const data = await apiUpload("/api/ora/genes/upload", formData);
                STATE.ora.geneTable = data.table;
                STATE.ora.uploadedFilename = file.name;
                STATE.ora.duplicatesRemoved = data.duplicates_removed;
                STATE.ora.resolved = null;
                STATE.ora.finalTable = null;
                STATE.ora.ambiguousChoices = {};
                STATE.ora.scoreColumn = "";
                STATE.ora.scoreSummary = null;
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));
    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Example file");
    exampleLink.addEventListener("click", () => {
        // Tailor the example to the currently selected identifier type, same
        // as the custom-database example does for its own input type.
        const fname = idType === "ensembl_gene_id" ? "example_burden_input_ensembl.csv" : "example_burden_input_symbol.csv";
        downloadFile("/api/burden/example?identifier_type=" + encodeURIComponent(idType), fname);
    });
    wrap.appendChild(el("div", { class: "caption" }, [
        "Required columns: " + idColLabel + " and a numeric Variant Score (e.g. a CADD score, " +
        "a GWAS p-value, or another per-gene severity/association metric you provide). " +
        "Optional column: UniProt ID. ",
        exampleLink,
    ]));

    if (STATE.ora.geneTable) {
        if (STATE.ora.duplicatesRemoved > 0) {
            wrap.appendChild(alertBox("info", "Removed " + STATE.ora.duplicatesRemoved + " duplicate gene(s); kept the most informative occurrence of each."));
        }
        const maxHeight = Math.min(400, Math.max(200, (STATE.ora.geneTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.ora.geneTable, { maxHeight }));
        wrap.appendChild(renderBurdenScoreColumnWidget());
    }

    return wrap;
}

// ---------------------------------------------------------------------------
// Gene Function Annotation: retrieval screen (table mode + single-gene mode).
// Table mode reuses the shared upload/resolve flow (like ORA - no score
// column needed). Single-gene mode writes STATE.ora.finalTable directly with
// one row, mirroring what the server-side resolve flow produces, so the
// existing renderContinueToAnalysisButton() gate works unchanged for it too.
// ---------------------------------------------------------------------------

function renderGeneFunctionAnnotationScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene Function Annotation"));

    const typeRow = el("div", { class: "row" });
    function typeBtn(label, type) {
        const btn = el("button", { class: "btn" }, label);
        if (STATE.ora.identifierType === type) btn.style.fontWeight = "800";
        btn.addEventListener("click", () => {
            if (STATE.ora.identifierType !== type) {
                STATE.ora = freshOraState();
                STATE.ora.identifierType = type;
                render();
            }
        });
        return btn;
    }
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Gene Symbol", "gene_symbol")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Ensembl Gene ID", "ensembl_gene_id")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Single gene", "single_gene")));
    typeRow.appendChild(el("div", { class: "col", style: { flex: "1" } }));
    wrap.appendChild(typeRow);

    if (!STATE.ora.identifierType) {
        wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
        wrap.appendChild(backToGenomicsActions());
        return wrap;
    }

    if (STATE.ora.identifierType === "single_gene") {
        wrap.appendChild(renderAnnotationSingleGeneSection());
    } else {
        if (STATE.ora.identifierType === "gene_symbol") {
            wrap.appendChild(el("div", { class: "section-label" }, "Select organism"));
            wrap.appendChild(renderOraOrganismSelector());
        }
        wrap.appendChild(renderAnnotationUploadSection());
        wrap.appendChild(renderGeneResolveControls());
        if (STATE.ora.resolved) {
            wrap.appendChild(renderGeneResults());
        }
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderAnnotationUploadSection() {
    const wrap = el("div", {});
    const idType = STATE.ora.identifierType;
    const idColLabel = idType === "gene_symbol" ? "Gene Symbol" : "Ensembl Gene ID";

    wrap.appendChild(el("div", { class: "section-label" }, "Upload gene table"));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.ora.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            formData.append("identifier_type", idType);
            try {
                const data = await apiUpload("/api/ora/genes/upload", formData);
                STATE.ora.geneTable = data.table;
                STATE.ora.uploadedFilename = file.name;
                STATE.ora.duplicatesRemoved = data.duplicates_removed;
                STATE.ora.resolved = null;
                STATE.ora.finalTable = null;
                STATE.ora.ambiguousChoices = {};
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));
    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Example file");
    exampleLink.addEventListener("click", () => {
        // Reuses the Burden Test's example files (same required columns; the
        // extra Variant Score column those carry is simply irrelevant here).
        const fname = idType === "ensembl_gene_id" ? "example_burden_input_ensembl.csv" : "example_burden_input_symbol.csv";
        downloadFile("/api/annotation/example?identifier_type=" + encodeURIComponent(idType), fname);
    });
    wrap.appendChild(el("div", { class: "caption" }, [
        "Required column: " + idColLabel + ". Optional column: UniProt ID. ",
        exampleLink,
    ]));

    if (STATE.ora.geneTable) {
        if (STATE.ora.duplicatesRemoved > 0) {
            wrap.appendChild(alertBox("info", "Removed " + STATE.ora.duplicatesRemoved + " duplicate gene(s); kept the most informative occurrence of each."));
        }
        const maxHeight = Math.min(400, Math.max(200, (STATE.ora.geneTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.ora.geneTable, { maxHeight }));
    }

    return wrap;
}

function renderAnnotationSingleGeneSection() {
    const wrap = el("div", {});
    wrap.appendChild(el("div", { class: "section-label" }, "Look up a single gene"));
    wrap.appendChild(el("div", { class: "caption" },
        "Enter a Gene Symbol, Ensembl Gene ID, or UniProt accession directly. A UniProt accession " +
        "(e.g. P38398) skips resolution entirely; a Gene Symbol needs an organism selected below " +
        "(an Ensembl Gene ID is species-specific and resolves without one)."));

    wrap.appendChild(el("div", { class: "section-label" }, "Organism (only needed for a Gene Symbol)"));
    wrap.appendChild(renderOraOrganismSelector());

    const row = el("div", { class: "row", style: { alignItems: "center" } });
    const input = el("input", { class: "text-input", type: "text", placeholder: "e.g. BRCA1, ENSG00000012048, or P38398" });
    input.value = STATE.ora.annotationSingleIdentifier || "";
    input.addEventListener("input", () => { STATE.ora.annotationSingleIdentifier = input.value; });
    row.appendChild(el("div", { class: "col" }, input));

    const btn = el("button", { class: "btn" }, "Look up");
    btn.addEventListener("click", async () => {
        const identifier = (STATE.ora.annotationSingleIdentifier || "").trim();
        if (!identifier) { alert("Enter a gene identifier."); return; }
        btn.disabled = true;
        try {
            const data = await apiPostJson("/api/annotation/single", { identifier, taxon_id: STATE.ora.taxonId || "" });
            STATE.ora.annotationSingleResult = data;
            STATE.ora.finalTable = {
                columns: ["Gene Symbol", "UniProt ID"],
                rows: [{ "Gene Symbol": identifier, "UniProt ID": data.uniprot_id }],
            };
            render();
        } catch (err) {
            STATE.ora.annotationSingleResult = null;
            STATE.ora.finalTable = null;
            alert(err.message);
        } finally {
            btn.disabled = false;
        }
    });
    row.appendChild(el("div", { class: "col", style: { flex: "0 0 auto" } }, btn));
    wrap.appendChild(row);

    if (STATE.ora.annotationSingleResult) {
        const res = STATE.ora.annotationSingleResult;
        const msg = res.skipped_resolution
            ? "Using UniProt accession " + res.uniprot_id + " directly (resolution skipped)."
            : "Resolved to UniProt accession " + res.uniprot_id + ".";
        wrap.appendChild(alertBox("info", msg));
        wrap.appendChild(renderContinueToAnalysisButton());
    }

    return wrap;
}

// Known aliases used to auto-detect the variant-score column so most users
// never have to touch the dropdown.
const BURDEN_SCORE_COLUMN_HINTS = ["variant score", "score", "cadd", "cadd score", "p-value", "pvalue", "p value", "severity"];

function renderBurdenScoreColumnWidget() {
    const idCols = ["gene symbol", "ensembl gene id", "uniprot id"];
    const candidateCols = (STATE.ora.geneTable.columns || []).filter(
        (c) => idCols.indexOf(c.trim().toLowerCase()) === -1
    );

    const wrap = el("div", { class: "ora-score-column-widget" });
    wrap.appendChild(el("div", { class: "section-label" }, "Variant score column"));
    wrap.appendChild(el("div", { class: "caption" },
        "The value is read as a plain number - it can be any metric (a CADD score, a raw GWAS p-value, etc.); " +
        "the tool doesn't need to know which. The Direction below tells it how to orient that number. " +
        "For a raw p-value, pick \"Lower = more burden\" directly - no need to convert it to -log10(p) yourself, " +
        "since the test only depends on rank order, which negation preserves identically."));

    if (candidateCols.length === 0) {
        wrap.appendChild(alertBox("warning", "No extra column found to use as the variant score. Add a numeric column (e.g. \"Variant Score\") to your file and re-upload."));
        return wrap;
    }

    if (!STATE.ora.scoreColumn || candidateCols.indexOf(STATE.ora.scoreColumn) === -1) {
        const hinted = candidateCols.find((c) => BURDEN_SCORE_COLUMN_HINTS.indexOf(c.trim().toLowerCase()) !== -1);
        STATE.ora.scoreColumn = hinted || candidateCols[0];
    }

    const row = el("div", { class: "ora-score-column-row" });
    row.appendChild(el("span", { class: "chart-option-label" }, "Column:"));
    const colSelect = el("select", { class: "chart-option-select" });
    candidateCols.forEach((c) => colSelect.appendChild(el("option", { value: c }, c)));
    colSelect.value = STATE.ora.scoreColumn;
    colSelect.addEventListener("change", () => { STATE.ora.scoreColumn = colSelect.value; STATE.ora.scoreSummary = null; });
    row.appendChild(colSelect);

    row.appendChild(el("span", { class: "chart-option-label" }, "Direction:"));
    const dirSelect = el("select", { class: "chart-option-select" });
    dirSelect.appendChild(el("option", { value: "higher" }, "Higher = more burden"));
    dirSelect.appendChild(el("option", { value: "lower" }, "Lower = more burden"));
    dirSelect.value = STATE.ora.scoreDirection || "higher";
    dirSelect.addEventListener("change", () => { STATE.ora.scoreDirection = dirSelect.value; STATE.ora.scoreSummary = null; });
    row.appendChild(dirSelect);

    const confirmBtn = el("button", { class: "btn btn-small", type: "button" }, "Set score column");
    confirmBtn.addEventListener("click", async () => {
        try {
            const data = await apiPostJson("/api/ora/genes/score-column", {
                column: STATE.ora.scoreColumn,
                direction: STATE.ora.scoreDirection || "higher",
            });
            STATE.ora.scoreSummary = {
                column: STATE.ora.scoreColumn,
                direction: STATE.ora.scoreDirection || "higher",
                n_scored: data.n_scored,
                n_non_numeric: data.n_non_numeric,
            };
            render();
        } catch (err) {
            STATE.ora.scoreSummary = null;
            alert(err.message);
            render();
        }
    });
    row.appendChild(confirmBtn);
    wrap.appendChild(row);

    if (STATE.ora.scoreSummary) {
        const s = STATE.ora.scoreSummary;
        let msg = "Using \"" + s.column + "\" (" + (s.direction === "lower" ? "lower = more burden" : "higher = more burden") + ") - " + s.n_scored + " gene(s) scored.";
        if (s.n_non_numeric > 0) msg += " " + s.n_non_numeric + " row(s) had a non-numeric value and were skipped.";
        wrap.appendChild(el("div", { class: "caption ora-score-summary" }, msg));
    } else {
        wrap.appendChild(el("div", { class: "caption" }, "Confirm the column and direction, then click \"Set score column\" before continuing to analysis."));
    }

    return wrap;
}

// ---------------------------------------------------------------------------
// Network Propagation / Module Detection: retrieval screen.
// Retrieval is shared with ORA/Burden (same STATE.ora, resolve flow), with two
// module-specific differences: (1) the organism selector is driven by STRING's
// OWN species list (STRING covers a subset of taxa - picking from STRING
// guarantees a network exists), and it is shown for BOTH identifier types
// because STRING needs the taxon regardless of how genes are resolved; (2) the
// score column is OPTIONAL (sets seed strength; omitted -> equal-weight seeds).
// ---------------------------------------------------------------------------

// STRING-species organism selector. Sets STATE.ora.taxonId/taxonName (the NCBI
// taxon, which UniProt gene resolution also understands) so one pick serves
// both STRING and resolution. Mirrors renderOraOrganismSelector but queries
// /api/netprop/species instead of the UniProt taxonomy search.
function renderNetpropOrganismSelector() {
    const wrap = el("div", {});
    const input = el("input", { class: "text-input", type: "text", placeholder: "Species / organism in STRING (name or NCBI taxon ID)" });
    const resultsSelect = el("select", { class: "hidden" });
    const selectedCaption = el("div", { class: "caption" });

    function refreshCaption() {
        selectedCaption.textContent = STATE.ora.taxonId ? "Selected: " + STATE.ora.taxonName + " (taxon " + STATE.ora.taxonId + ")" : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = await apiJson("/api/netprop/species?query=" + encodeURIComponent(query));
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                (data.results || []).forEach((rec, idx) => {
                    resultsSelect.appendChild(el("option", { value: String(idx) }, rec.taxon_id + " | " + rec.name));
                });
                resultsSelect._records = data.results || [];
                resultsSelect.classList.toggle("hidden", (data.results || []).length === 0);
            } catch (err) { /* keep prior dropdown on transient errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        STATE.ora.taxonId = String(rec.taxon_id || "");
        STATE.ora.taxonName = String(rec.name || "");
        render();
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function renderGenomicsNetpropScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-module-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Network Propagation / Module Detection"));

    const typeRow = el("div", { class: "row" });
    function typeBtn(label, type) {
        const btn = el("button", { class: "btn" }, label);
        if (STATE.ora.identifierType === type) btn.style.fontWeight = "800";
        btn.addEventListener("click", () => {
            if (STATE.ora.identifierType !== type) {
                STATE.ora = freshOraState();
                STATE.ora.identifierType = type;
                render();
            }
        });
        return btn;
    }
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Gene Symbol", "gene_symbol")));
    typeRow.appendChild(el("div", { class: "col" }, typeBtn("Ensembl Gene ID", "ensembl_gene_id")));
    typeRow.appendChild(el("div", { class: "col", style: { flex: "1" } }));
    wrap.appendChild(typeRow);

    if (!STATE.ora.identifierType) {
        wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
        wrap.appendChild(backToGenomicsActions());
        return wrap;
    }

    // Organism selector shown for BOTH identifier types - STRING needs the
    // taxon even when Ensembl resolution itself doesn't.
    wrap.appendChild(el("div", { class: "section-label" }, "Select organism (STRING network)"));
    wrap.appendChild(renderNetpropOrganismSelector());

    wrap.appendChild(renderNetpropUploadSection());
    wrap.appendChild(renderGeneResolveControls());

    if (STATE.ora.resolved) {
        wrap.appendChild(renderGeneResults());
    }

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderNetpropUploadSection() {
    const wrap = el("div", {});
    const idType = STATE.ora.identifierType;
    const idColLabel = idType === "gene_symbol" ? "Gene Symbol" : "Ensembl Gene ID";

    wrap.appendChild(el("div", { class: "section-label" }, "Upload gene table"));
    wrap.appendChild(fileInputWidget({
        accept: ".csv",
        buttonLabel: "Choose CSV file",
        initialFilename: STATE.ora.uploadedFilename,
        onChange: async (file) => {
            if (!file) return;
            const formData = new FormData();
            formData.append("file", file);
            formData.append("identifier_type", idType);
            try {
                const data = await apiUpload("/api/ora/genes/upload", formData);
                STATE.ora.geneTable = data.table;
                STATE.ora.uploadedFilename = file.name;
                STATE.ora.duplicatesRemoved = data.duplicates_removed;
                STATE.ora.resolved = null;
                STATE.ora.finalTable = null;
                STATE.ora.ambiguousChoices = {};
                STATE.ora.scoreColumn = "";
                STATE.ora.scoreSummary = null;
                render();
            } catch (err) {
                alert(err.message);
            }
        },
    }));
    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "Example file");
    exampleLink.addEventListener("click", () => {
        const fname = idType === "ensembl_gene_id" ? "example_netprop_input_ensembl.csv" : "example_netprop_input_symbol.csv";
        downloadFile("/api/netprop/example?identifier_type=" + encodeURIComponent(idType), fname);
    });
    wrap.appendChild(el("div", { class: "caption" }, [
        "Required column: " + idColLabel + ". Optional columns: UniProt ID, and a numeric Gene Score " +
        "(e.g. a CADD score or GWAS metric) that sets each seed's starting strength - omit it to treat all " +
        "genes as equal-weight seeds. ",
        exampleLink,
    ]));

    if (STATE.ora.geneTable) {
        if (STATE.ora.duplicatesRemoved > 0) {
            wrap.appendChild(alertBox("info", "Removed " + STATE.ora.duplicatesRemoved + " duplicate gene(s); kept the most informative occurrence of each."));
        }
        const maxHeight = Math.min(400, Math.max(200, (STATE.ora.geneTable.rows.length + 1) * 35));
        wrap.appendChild(buildTable(STATE.ora.geneTable, { maxHeight }));
        // Optional score-column widget: only when the file has an extra numeric
        // column to offer as the seed weight. Reuses the Burden widget, which
        // stores STATE.ora.scoreSummary (consumed by _job_run_netprop as seed
        // strengths). Skipping it is fine (equal-weight seeds).
        const idCols = ["gene symbol", "ensembl gene id", "uniprot id"];
        const hasExtra = (STATE.ora.geneTable.columns || []).some(c => idCols.indexOf(c.trim().toLowerCase()) === -1);
        if (hasExtra) {
            wrap.appendChild(el("div", { class: "caption ora-score-optional-note" },
                "Optional: pick a Gene Score column below to weight seeds; otherwise all genes seed equally."));
            wrap.appendChild(renderBurdenScoreColumnWidget());
        }
    }

    return wrap;
}

// ---------------------------------------------------------------------------
// ORA statistical analysis: configure -> run -> results
// ---------------------------------------------------------------------------

function renderOraAnalysisConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Over-Representation Analysis: database"));

    const locked = STATE.ora.analysisRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const panel = el("div", { class: "ora-config-panel" });

    const availableDbs = ORA_AVAILABLE_DATABASES.filter(
        db => db !== "GWAS Catalog" || STATE.ora.identifierType === "gene_symbol"
    );
    availableDbs.forEach((dbName) => {
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked;
        cb.checked = STATE.ora.analysisConfig.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (STATE.ora.analysisConfig.databases.indexOf(dbName) === -1) {
                    STATE.ora.analysisConfig.databases.push(dbName);
                }
                if (!STATE.ora.analysisConfig.perDatabase[dbName]) {
                    STATE.ora.analysisConfig.perDatabase[dbName] = { taxonId: "", taxonName: "" };
                }
            } else {
                STATE.ora.analysisConfig.databases = STATE.ora.analysisConfig.databases.filter((d) => d !== dbName);
            }
            renderOraConfigPanel(panel);
            renderOraRunButton(runBtnWrap);
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    wrap.appendChild(panel);
    renderOraConfigPanel(panel);

    const actions = el("div", { class: "right-actions ora-config-actions" });
    const btnRow = el("div", { class: "btn-row ora-config-btn-row" });
    const backBtn = el("button", { class: "btn ora-back-btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "genomics_module"; render(); });
    btnRow.appendChild(backBtn);
    const runBtnWrap = el("div", { class: "col" });
    renderOraRunButton(runBtnWrap);
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function _oraConfigIsRunnable() {
    const cfg = STATE.ora.analysisConfig;
    if (cfg.databases.length === 0) return false;
    return cfg.databases.every((d) => {
        const c = cfg.perDatabase[d];
        if (!c) return false;
        if (d === ORA_CUSTOM_DB_LABEL) return !!c.uploaded;
        return !!c.taxonId;
    });
}

function renderOraRunButton(container) {
    container.innerHTML = "";
    const locked = STATE.ora.analysisRunLocked;
    const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
    runBtn.disabled = locked || !_oraConfigIsRunnable();
    runBtn.addEventListener("click", startOraAnalysisRun);
    container.appendChild(runBtn);
}

const ORA_DEFAULT_SIG_THRESHOLD = 0.05;

// Shared "Significance threshold (FDR <)" control, added to every database block
// (including Custom database). Stored per-database so each result block/chart can
// use its own threshold; display-only client-side, never sent to the server.
function renderOraSigThresholdRow(dbName, locked, configRoot) {
    configRoot = configRoot || STATE.ora.analysisConfig;
    const dbCfg = configRoot.perDatabase[dbName] || {};
    configRoot.perDatabase[dbName] = dbCfg;
    if (!dbCfg.sigThreshold && dbCfg.sigThreshold !== 0) dbCfg.sigThreshold = ORA_DEFAULT_SIG_THRESHOLD;

    const row = el("div", { class: "ora-sig-threshold-row" });
    row.appendChild(el("span", { class: "chart-option-label" }, "Significance threshold (FDR <)"));
    const input = el("input", {
        type: "number", class: "ora-sig-threshold-input",
        min: "0", max: "1", step: "0.001",
    });
    input.value = dbCfg.sigThreshold;
    input.disabled = locked;
    input.addEventListener("change", () => {
        let v = parseFloat(input.value);
        if (!isFinite(v) || v <= 0 || v > 1) v = ORA_DEFAULT_SIG_THRESHOLD;
        dbCfg.sigThreshold = v;
        input.value = v;
    });
    row.appendChild(input);
    return row;
}

const ORA_DEFAULT_MIN_TERM_SIZE = 3;
// clusterProfiler's own maxGSSize default. GO ancestor propagation means
// near-root terms (e.g. "biological_process" itself) can now carry tens of
// thousands of genes; left uncapped, those would report as FDR-significant,
// fold~1.0 "enrichment" - statistically meaningless noise crowding out real
// hits. Still clearable to blank for "no limit" if a run genuinely needs it.
const ORA_DEFAULT_MAX_TERM_SIZE = 500;

// Optional term-size filter (all databases). min (default 3) excludes tiny,
// statistically weak sets; max (default 500, blank = no limit) excludes huge
// generic terms (e.g. "protein binding", or a propagated near-root term).
// Stored per-database and sent to the server, which re-runs Fisher/BH over
// only the surviving terms.
function renderOraTermSizeRow(dbName, locked, configRoot) {
    configRoot = configRoot || STATE.ora.analysisConfig;
    const dbCfg = configRoot.perDatabase[dbName] || {};
    configRoot.perDatabase[dbName] = dbCfg;
    if (!dbCfg.minTermSize && dbCfg.minTermSize !== 0) dbCfg.minTermSize = ORA_DEFAULT_MIN_TERM_SIZE;
    if (typeof dbCfg.maxTermSize === "undefined") dbCfg.maxTermSize = ORA_DEFAULT_MAX_TERM_SIZE;

    const row = el("div", { class: "ora-term-size-row" });
    row.appendChild(el("span", { class: "chart-option-label" }, "Term size (genes)"));

    const minInput = el("input", {
        type: "number", class: "ora-term-size-input", min: "1", step: "1",
        placeholder: "min", title: "Minimum genes a term must have to be tested",
    });
    minInput.value = dbCfg.minTermSize;
    minInput.disabled = locked;
    minInput.addEventListener("change", () => {
        let v = parseInt(minInput.value, 10);
        if (!isFinite(v) || v < 1) v = ORA_DEFAULT_MIN_TERM_SIZE;
        dbCfg.minTermSize = v;
        minInput.value = v;
    });

    const maxInput = el("input", {
        type: "number", class: "ora-term-size-input", min: "1", step: "1",
        placeholder: "max (" + ORA_DEFAULT_MAX_TERM_SIZE + ")", title: "Maximum genes a term may have (default " + ORA_DEFAULT_MAX_TERM_SIZE + "; blank = no limit)",
    });
    maxInput.value = dbCfg.maxTermSize;
    maxInput.disabled = locked;
    maxInput.addEventListener("change", () => {
        const raw = maxInput.value.trim();
        if (raw === "") { dbCfg.maxTermSize = ""; return; }
        let v = parseInt(raw, 10);
        if (!isFinite(v) || v < 1) { dbCfg.maxTermSize = ""; maxInput.value = ""; return; }
        dbCfg.maxTermSize = v;
        maxInput.value = v;
    });

    row.appendChild(minInput);
    row.appendChild(el("span", { class: "ora-term-size-dash" }, "-"));
    row.appendChild(maxInput);
    return row;
}

// Custom background (reference) gene list - GO/Reactome only (accession space).
// Restricts the ORA universe to |annotated universe intersect background|. Optional;
// leaving it empty (checkbox off) uses the whole annotated genome. Gated behind
// a checkbox (same pattern as the IEA toggle) so the upload widget/explanation/
// summary stay hidden until the user opts in - keeps the default card compact.
function renderOraBackgroundWidget(dbName, locked) {
    const dbCfg = STATE.ora.analysisConfig.perDatabase[dbName] || {};
    STATE.ora.analysisConfig.perDatabase[dbName] = dbCfg;
    if (typeof dbCfg.useBackground === "undefined") dbCfg.useBackground = false;

    const wrap = el("div", { class: "ora-background-widget" });

    const toggleRow = el("div", { class: "ora-iea-row ora-background-toggle-row" });
    const toggleCb = el("input", { type: "checkbox", class: "ora-iea-checkbox ora-background-checkbox" });
    toggleCb.checked = dbCfg.useBackground;
    toggleCb.disabled = locked;
    toggleCb.addEventListener("change", () => {
        dbCfg.useBackground = toggleCb.checked;
        render();
    });
    toggleRow.appendChild(toggleCb);
    toggleRow.appendChild(document.createTextNode(" Background gene list (optional)"));
    wrap.appendChild(toggleRow);

    // Collapsed: uploaded background (if any) stays in memory but is inert - // has_background is only sent true when this checkbox is also on.
    if (!dbCfg.useBackground) return wrap;

    const body = el("div", { class: "ora-background-body" });

    if (!locked) {
        body.appendChild(fileInputWidget({
            accept: ".txt,.csv,.tsv",
            buttonLabel: "Choose background file",
            initialFilename: dbCfg.backgroundSummary ? dbCfg.backgroundSummary.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                formData.append("database", dbName);
                try {
                    const data = await apiUpload("/api/ora/analysis/background-upload", formData);
                    dbCfg.backgroundSummary = { filename: data.filename, n_background: data.n_background };
                    render();
                } catch (err) {
                    dbCfg.backgroundSummary = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    if (dbCfg.backgroundSummary) {
        const s = dbCfg.backgroundSummary;
        body.appendChild(el("div", { class: "caption ora-background-summary" },
            "Background: \"" + s.filename + "\" - " + s.n_background + " UniProt accessions. " +
            "The universe is restricted to genes present in both this list and the annotated database."));
        if (!locked) {
            const removeBtn = el("button", { class: "btn btn-small ora-background-remove", type: "button" }, "Remove background");
            removeBtn.addEventListener("click", () => {
                dbCfg.backgroundSummary = null;
                render();
            });
            body.appendChild(removeBtn);
        }
    } else {
        const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "example file");
        exampleLink.addEventListener("click", () => {
            downloadFile("/api/ora/background-example", "example_ora_background.csv");
        });
        body.appendChild(el("div", { class: "caption ora-background-help" }, [
            "Provide UniProt accessions (one per line, or a column headed \"UniProt ID\") to use a custom reference set " +
            "(e.g. genes expressed in your RNA-seq assay). Leave empty to use the whole annotated genome. (",
            exampleLink,
            ")",
        ]));
    }

    wrap.appendChild(body);
    return wrap;
}

function renderOraConfigPanel(panel) {
    panel.innerHTML = "";
    const locked = STATE.ora.analysisRunLocked;
    STATE.ora.analysisConfig.databases.forEach((dbName) => {
        const dbClass = dbName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
        const block = el("div", { class: "rank-card ora-db-card ora-db-card-" + dbClass });
        block.appendChild(el("div", { class: "section-label ora-db-title" }, dbName));
        if (ORA_DATABASE_CAPTIONS[dbName]) block.appendChild(el("div", { class: "caption ora-db-caption" }, ORA_DATABASE_CAPTIONS[dbName]));

        // Custom database uses a file-upload widget instead of a species selector.
        if (dbName === ORA_CUSTOM_DB_LABEL) {
            block.appendChild(renderOraCustomDbWidget(locked, panel));
            block.appendChild(renderOraTermSizeRow(dbName, locked));
            block.appendChild(renderOraSigThresholdRow(dbName, locked));
            panel.appendChild(block);
            return;
        }

        block.appendChild(renderOraSpeciesSearchWidget(dbName, locked, panel));

        // IEA toggle - only for GO and Reactome (GWAS has no evidence-code concept)
        if (dbName === "GO" || dbName === "Reactome") {
            const dbCfg = STATE.ora.analysisConfig.perDatabase[dbName] || {};
            if (typeof dbCfg.includeIea === "undefined") dbCfg.includeIea = true;
            const ieaRow = el("div", { class: "ora-iea-row" });
            const ieaCb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
            ieaCb.checked = dbCfg.includeIea;
            ieaCb.disabled = locked;
            ieaCb.addEventListener("change", () => {
                dbCfg.includeIea = ieaCb.checked;
            });
            ieaRow.appendChild(ieaCb);
            ieaRow.appendChild(document.createTextNode(" Include electronic (IEA) annotations"));
            block.appendChild(ieaRow);
        }

        // Custom background upload - GO/Reactome only (accession space).
        if (dbName === "GO" || dbName === "Reactome") {
            block.appendChild(renderOraBackgroundWidget(dbName, locked));
        }

        block.appendChild(renderOraTermSizeRow(dbName, locked));
        block.appendChild(renderOraSigThresholdRow(dbName, locked));

        panel.appendChild(block);
    });
}

// ---------------------------------------------------------------------------
// Gene-set Burden Test: configure -> run -> results. Reuses every ORA config
// widget above (species search, custom-DB widget, term-size row, sig-threshold
// row) via their configRoot parameter pointed at STATE.ora.burdenConfig instead
// of STATE.ora.analysisConfig. No background widget (not applicable to a
// competitive rank test) and no GWAS Catalog (see BURDEN_AVAILABLE_DATABASES).
// ---------------------------------------------------------------------------

function renderBurdenConfigPanel(panel) {
    panel.innerHTML = "";
    const locked = STATE.ora.burdenRunLocked;
    const configRoot = STATE.ora.burdenConfig;
    configRoot.databases.forEach((dbName) => {
        const dbClass = dbName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
        const block = el("div", { class: "rank-card ora-db-card ora-db-card-" + dbClass });
        block.appendChild(el("div", { class: "section-label ora-db-title" }, dbName));
        if (BURDEN_DATABASE_CAPTIONS[dbName]) block.appendChild(el("div", { class: "caption ora-db-caption" }, BURDEN_DATABASE_CAPTIONS[dbName]));

        if (dbName === ORA_CUSTOM_DB_LABEL) {
            block.appendChild(renderOraCustomDbWidget(locked, panel, configRoot));
            block.appendChild(renderOraTermSizeRow(dbName, locked, configRoot));
            block.appendChild(renderOraSigThresholdRow(dbName, locked, configRoot));
            panel.appendChild(block);
            return;
        }

        block.appendChild(renderOraSpeciesSearchWidget(dbName, locked, panel, configRoot));

        // IEA toggle - only for GO and Reactome (GWAS has no evidence-code concept)
        if (dbName === "GO" || dbName === "Reactome") {
            const dbCfg = configRoot.perDatabase[dbName] || {};
            if (typeof dbCfg.includeIea === "undefined") dbCfg.includeIea = true;
            const ieaRow = el("div", { class: "ora-iea-row" });
            const ieaCb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
            ieaCb.checked = dbCfg.includeIea;
            ieaCb.disabled = locked;
            ieaCb.addEventListener("change", () => {
                dbCfg.includeIea = ieaCb.checked;
            });
            ieaRow.appendChild(ieaCb);
            ieaRow.appendChild(document.createTextNode(" Include electronic (IEA) annotations"));
            block.appendChild(ieaRow);
        }

        block.appendChild(renderOraTermSizeRow(dbName, locked, configRoot));
        block.appendChild(renderOraSigThresholdRow(dbName, locked, configRoot));

        panel.appendChild(block);
    });
}

function renderBurdenAnalysisConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene-set Burden Test: database"));

    const locked = STATE.ora.burdenRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Configuration is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    wrap.appendChild(el("div", { class: "section-label ora-select-label" }, "Select database(s)"));
    const checkRow = el("div", { class: "db-checkbox-row ora-db-checkbox-row" });
    const panel = el("div", { class: "ora-config-panel" });

    BURDEN_AVAILABLE_DATABASES.forEach((dbName) => {
        const label = el("label", { class: "ora-db-pill" });
        const cb = el("input", { type: "checkbox", class: "ora-db-checkbox" });
        cb.disabled = locked;
        cb.checked = STATE.ora.burdenConfig.databases.indexOf(dbName) !== -1;
        cb.addEventListener("change", () => {
            if (cb.checked) {
                if (STATE.ora.burdenConfig.databases.indexOf(dbName) === -1) {
                    STATE.ora.burdenConfig.databases.push(dbName);
                }
                if (!STATE.ora.burdenConfig.perDatabase[dbName]) {
                    STATE.ora.burdenConfig.perDatabase[dbName] = { taxonId: "", taxonName: "" };
                }
            } else {
                STATE.ora.burdenConfig.databases = STATE.ora.burdenConfig.databases.filter((d) => d !== dbName);
            }
            renderBurdenConfigPanel(panel);
            renderBurdenRunButton(runBtnWrap);
        });
        label.appendChild(cb);
        label.appendChild(document.createTextNode(" " + dbName));
        checkRow.appendChild(label);
    });
    wrap.appendChild(checkRow);
    wrap.appendChild(panel);
    renderBurdenConfigPanel(panel);

    const actions = el("div", { class: "right-actions ora-config-actions" });
    const btnRow = el("div", { class: "btn-row ora-config-btn-row" });
    const backBtn = el("button", { class: "btn ora-back-btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "genomics_module"; render(); });
    btnRow.appendChild(backBtn);
    const runBtnWrap = el("div", { class: "col" });
    renderBurdenRunButton(runBtnWrap);
    btnRow.appendChild(runBtnWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    return wrap;
}

function _burdenConfigIsRunnable() {
    const cfg = STATE.ora.burdenConfig;
    if (cfg.databases.length === 0) return false;
    return cfg.databases.every((d) => {
        const c = cfg.perDatabase[d];
        if (!c) return false;
        if (d === ORA_CUSTOM_DB_LABEL) return !!c.uploaded;
        return !!c.taxonId;
    });
}

function renderBurdenRunButton(container) {
    container.innerHTML = "";
    const locked = STATE.ora.burdenRunLocked;
    const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
    runBtn.disabled = locked || !_burdenConfigIsRunnable();
    runBtn.addEventListener("click", startBurdenAnalysisRun);
    container.appendChild(runBtn);
}

function startBurdenAnalysisRun() {
    if (!_burdenConfigIsRunnable()) return;
    const perDbApi = {};
    STATE.ora.burdenConfig.databases.forEach(db => {
        const cfg = STATE.ora.burdenConfig.perDatabase[db] || {};
        const minTs = (cfg.minTermSize || cfg.minTermSize === 0) ? cfg.minTermSize : ORA_DEFAULT_MIN_TERM_SIZE;
        // "" is an explicit user opt-out (no cap); undefined (row never
        // rendered/touched) falls back to the same default the term-size
        // row itself uses, not to "no cap".
        const maxTs = cfg.maxTermSize === "" ? null : (typeof cfg.maxTermSize === "undefined" ? ORA_DEFAULT_MAX_TERM_SIZE : cfg.maxTermSize);
        if (db === ORA_CUSTOM_DB_LABEL) {
            perDbApi[db] = { id_type: cfg.idType || "accession", min_term_size: minTs, max_term_size: maxTs };
        } else {
            perDbApi[db] = {
                taxon_id: cfg.taxonId,
                taxon_name: cfg.taxonName,
                include_iea: typeof cfg.includeIea === "undefined" ? true : cfg.includeIea,
                min_term_size: minTs,
                max_term_size: maxTs,
            };
        }
    });
    const payload = {
        databases: STATE.ora.burdenConfig.databases,
        per_database: perDbApi,
    };
    apiPostJson("/api/burden/analysis/configure", payload).then(() => {
        STATE.ora.burdenRunLocked = true;
        STATE.ora.burdenResults = null;
        STATE.aaPage = "genomics_burden_results";
        render();
        return apiPostJson("/api/burden/analysis/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("ora-analysis-progress-inner");
                const label = document.getElementById("ora-analysis-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => {
                STATE.ora.burdenResults = result;
                render();
            },
            onError: (err) => {
                alert("Burden test run failed: " + err);
                STATE.ora.burdenRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.ora.burdenRunLocked = false;
        render();
    });
}

// ---------------------------------------------------------------------------
// Network Propagation / Module Detection: configure -> run -> results.
// The organism (STRING taxon) is already chosen on the retrieval screen; the
// config screen only sets the two checkboxes (module detection auto-implies
// propagation) and the exposed parameters.
// ---------------------------------------------------------------------------

const NETPROP_CONFIDENCE_CHOICES = [["0.15 (low)", "150"], ["0.4 (medium)", "400"], ["0.7 (high)", "700"], ["0.9 (highest)", "900"]];
const NETPROP_PERMUTATION_CHOICES = [["500", "500"], ["1000", "1000"], ["2000", "2000"], ["5000", "5000"], ["10000", "10000"]];

function renderNetpropAnalysisConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Network Propagation / Module Detection: setup"));

    const locked = STATE.ora.netpropRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Setup is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    const cfg = STATE.ora.netpropConfig;
    wrap.appendChild(el("div", { class: "caption" },
        "Organism: " + (STATE.ora.taxonName || "?") + " (STRING taxon " + (STATE.ora.taxonId || "?") + "). " +
        (STATE.ora.scoreSummary ? "Seeds weighted by \"" + STATE.ora.scoreSummary.column + "\"." : "Seeds are equal-weight (no score column).")));

    // --- Analysis stages (checkboxes) ---
    wrap.appendChild(el("div", { class: "section-label" }, "Analysis"));
    const stageBox = el("div", { class: "rank-card ora-db-card" });
    const propRow = el("div", { class: "ora-iea-row" });
    const propCb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
    propCb.checked = true; propCb.disabled = true;   // always runs; it's the base analysis
    propRow.appendChild(propCb);
    propRow.appendChild(document.createTextNode(" Network propagation (always runs - spreads seed signal, keeps significant hits)"));
    stageBox.appendChild(propRow);

    const modRow = el("div", { class: "ora-iea-row" });
    const modCb = el("input", { type: "checkbox", class: "ora-iea-checkbox netprop-module-checkbox" });
    modCb.checked = !!cfg.moduleDetection; modCb.disabled = locked;
    modCb.addEventListener("change", () => { cfg.moduleDetection = modCb.checked; });
    modRow.appendChild(modCb);
    modRow.appendChild(document.createTextNode(" Module detection (cluster the hits into modules + label each by enrichment - auto-runs propagation first)"));
    stageBox.appendChild(modRow);
    wrap.appendChild(stageBox);

    // --- Parameters ---
    wrap.appendChild(el("div", { class: "section-label" }, "Parameters"));
    const paramBox = el("div", { class: "rank-card ora-db-card" });

    function selectRow(label, choices, current, onChange) {
        const row = el("div", { class: "ora-term-size-row" });
        row.appendChild(el("span", { class: "chart-option-label" }, label));
        const sel = el("select", { class: "chart-option-select" });
        choices.forEach(([lbl, val]) => sel.appendChild(el("option", { value: val }, lbl)));
        sel.value = String(current);
        sel.disabled = locked;
        sel.addEventListener("change", () => onChange(sel.value));
        row.appendChild(sel);
        return row;
    }
    function numberRow(label, value, step, min, max, onChange) {
        const row = el("div", { class: "ora-term-size-row" });
        row.appendChild(el("span", { class: "chart-option-label" }, label));
        const inp = el("input", { type: "number", class: "ora-term-size-input", step: String(step), min: String(min) });
        inp.value = value; inp.disabled = locked;
        inp.addEventListener("change", () => onChange(inp.value));
        row.appendChild(inp);
        return row;
    }

    paramBox.appendChild(selectRow("STRING confidence cutoff", NETPROP_CONFIDENCE_CHOICES, cfg.confidence, (v) => { cfg.confidence = Number(v); }));
    paramBox.appendChild(numberRow("Significance threshold (FDR <)", cfg.fdrAlpha, 0.001, 0, (v) => {
        let x = parseFloat(v); cfg.fdrAlpha = (isFinite(x) && x > 0 && x <= 1) ? x : 0.05;
    }));
    paramBox.appendChild(numberRow("Minimum module size", cfg.minModuleSize, 1, 2, (v) => {
        let x = parseInt(v, 10); cfg.minModuleSize = (isFinite(x) && x >= 2) ? x : 3;
    }));
    paramBox.appendChild(selectRow("Number of permutations", NETPROP_PERMUTATION_CHOICES, cfg.nPermutations, (v) => { cfg.nPermutations = Number(v); }));
    paramBox.appendChild(el("div", { class: "caption ora-custom-help" },
        "Fixed (reported with results for reproducibility): restart probability 0.5, MCL inflation 2.0, degree-preserving permutations, STRING v12.0."));
    wrap.appendChild(paramBox);

    const actions = el("div", { class: "right-actions ora-config-actions" });
    const btnRow = el("div", { class: "btn-row ora-config-btn-row" });
    const backBtn = el("button", { class: "btn ora-back-btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "genomics_module"; render(); });
    btnRow.appendChild(backBtn);
    const runWrap = el("div", { class: "col" });
    const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
    runBtn.disabled = locked || !STATE.ora.taxonId;
    runBtn.addEventListener("click", startNetpropAnalysisRun);
    runWrap.appendChild(runBtn);
    btnRow.appendChild(runWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function startNetpropAnalysisRun() {
    if (!STATE.ora.taxonId) return;
    const cfg = STATE.ora.netpropConfig;
    const payload = {
        taxon_id: STATE.ora.taxonId,
        taxon_name: STATE.ora.taxonName,
        module_detection: !!cfg.moduleDetection,
        confidence: cfg.confidence,
        fdr_alpha: cfg.fdrAlpha,
        min_module_size: cfg.minModuleSize,
        n_permutations: cfg.nPermutations,
    };
    apiPostJson("/api/netprop/analysis/configure", payload).then(() => {
        STATE.ora.netpropRunLocked = true;
        STATE.ora.netpropResults = null;
        STATE.aaPage = "genomics_netprop_results";
        render();
        return apiPostJson("/api/netprop/analysis/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("ora-analysis-progress-inner");
                const label = document.getElementById("ora-analysis-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => { STATE.ora.netpropResults = result; render(); },
            onError: (err) => {
                alert("Network propagation run failed: " + err);
                STATE.ora.netpropRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.ora.netpropRunLocked = false;
        render();
    });
}

function renderOraCustomDbWidget(locked, panel, configRoot) {
    configRoot = configRoot || STATE.ora.analysisConfig;
    const wrap = el("div", { class: "ora-custom-db-widget" });
    const dbCfg = configRoot.perDatabase[ORA_CUSTOM_DB_LABEL] || { idType: "accession", uploaded: false, summary: null };
    if (!dbCfg.idType) dbCfg.idType = "accession";
    configRoot.perDatabase[ORA_CUSTOM_DB_LABEL] = dbCfg;

    // ID-type selector - restricted to the identifier space(s) the current input
    // type can actually match: gene symbols only exist for Gene Symbol input,
    // Ensembl IDs only for Ensembl input, and BLAST/sequence resolves to UniProt
    // accessions only. Accession is always available.
    const idOptionsByInput = {
        gene_symbol: [["UniProt accession", "accession"], ["Gene symbol", "symbol"]],
        ensembl_gene_id: [["UniProt accession", "accession"], ["Ensembl Gene ID", "ensembl"]],
        sequence: [["UniProt accession", "accession"]],
    };
    const idOptions = idOptionsByInput[STATE.ora.identifierType] || [["UniProt accession", "accession"]];
    // Guard the stored value against options the current input type doesn't support.
    if (!idOptions.some(([, val]) => val === dbCfg.idType)) {
        dbCfg.idType = "accession";
        dbCfg.uploaded = false;
        dbCfg.summary = null;
    }

    const idRow = el("label", { class: "ora-custom-idtype-row" });
    idRow.appendChild(el("span", { class: "chart-option-label" }, "Gene identifiers in the file:"));
    const idSelect = el("select", { class: "chart-option-select" });
    idOptions.forEach(([lbl, val]) => {
        idSelect.appendChild(el("option", { value: val }, lbl));
    });
    idSelect.value = dbCfg.idType;
    idSelect.disabled = locked || idOptions.length < 2;
    idSelect.addEventListener("change", () => {
        dbCfg.idType = idSelect.value;
        // id_type is baked into the parsed file server-side; force a re-upload.
        if (dbCfg.uploaded) { dbCfg.uploaded = false; dbCfg.summary = null; render(); }
    });
    idRow.appendChild(idSelect);
    wrap.appendChild(idRow);

    if (STATE.ora.identifierType === "sequence") {
        wrap.appendChild(el("div", { class: "caption ora-custom-help" },
            "Homology (BLAST) input resolves your sequences to UniProt accessions, so a custom database must use UniProt accessions."));
    }

    if (!locked) {
        wrap.appendChild(fileInputWidget({
            accept: ".gmt,.txt,.tsv,.csv",
            buttonLabel: "Choose gene-set file",
            initialFilename: dbCfg.summary ? dbCfg.summary.filename : "",
            onChange: async (file) => {
                if (!file) return;
                const formData = new FormData();
                formData.append("file", file);
                formData.append("id_type", dbCfg.idType);
                try {
                    const data = await apiUpload("/api/ora/analysis/custom-db-upload", formData);
                    dbCfg.summary = data;
                    dbCfg.uploaded = true;
                    render();
                } catch (err) {
                    dbCfg.uploaded = false;
                    dbCfg.summary = null;
                    alert(err.message);
                    render();
                }
            },
        }));
    }

    const ID_TYPE_LABELS = { accession: "UniProt accessions", symbol: "gene symbols", ensembl: "Ensembl Gene IDs" };

    if (dbCfg.summary) {
        const s = dbCfg.summary;
        wrap.appendChild(el("div", { class: "caption ora-custom-summary" },
            "Loaded \"" + s.filename + "\": " + s.n_terms + " gene set(s), " + s.universe_size +
            " unique genes (" + (ID_TYPE_LABELS[s.id_type] || s.id_type) + ")."));
    }

    const exampleLink = el("button", { class: "ora-inline-example-link", type: "button" }, "(example input files)");
    exampleLink.addEventListener("click", () => {
        // Tailor the example to the ORA input type (what determines the valid
        // identifier space), so ensembl/homology don't get the gene-symbol example.
        downloadFile("/api/ora/custom-db-example?input_type=" + encodeURIComponent(STATE.ora.identifierType || "gene_symbol"), "custom_geneset_db_example.zip");
    });
    wrap.appendChild(el("div", { class: "caption ora-custom-help" }, [
        "Format: GMT (Term | Description | Gene1 ...) or a tidy table with columns term_id, term_name, gene. " +
        "Genes must match the identifier type selected above; the background universe is the union of all genes in the file. ",
        exampleLink,
    ]));

    return wrap;
}

function renderOraSpeciesSearchWidget(dbName, locked, panel, configRoot) {
    configRoot = configRoot || STATE.ora.analysisConfig;
    const wrap = el("div", { class: "ora-species-widget" });
    const dbCfg = configRoot.perDatabase[dbName] || { taxonId: "", taxonName: "" };
    configRoot.perDatabase[dbName] = dbCfg;

    const input = el("input", { class: "text-input ora-species-input", type: "text", placeholder: dbName + " species (name or UniProt taxon ID)" });
    input.disabled = locked;
    const resultsSelect = el("select", { class: "hidden ora-species-select" });
    const selectedCaption = el("div", { class: "caption ora-selected-caption" });

    function refreshCaption() {
        selectedCaption.textContent = dbCfg.taxonId ? "Selected: " + dbCfg.taxonName + " (taxon " + dbCfg.taxonId + ")" : "";
    }
    refreshCaption();

    let lastQuery = "";
    let timer = null;
    input.addEventListener("input", () => {
        const query = input.value.trim();
        if (timer) clearTimeout(timer);
        if (query.length < 3) {
            resultsSelect.innerHTML = "";
            resultsSelect.classList.add("hidden");
            return;
        }
        timer = setTimeout(async () => {
            if (query === lastQuery) return;
            lastQuery = query;
            try {
                const data = await apiJson("/api/contexts/search?database=" + encodeURIComponent("Rhea + UniProt") + "&query=" + encodeURIComponent(query) + "&limit=10");
                resultsSelect.innerHTML = "";
                resultsSelect.appendChild(el("option", { value: "" }, "Select organism"));
                (data.results || []).forEach((rec, idx) => {
                    const label = rec.context_id && rec.context_name ? rec.context_id + " | " + rec.context_name : (rec.context_name || rec.context_id || "");
                    resultsSelect.appendChild(el("option", { value: String(idx) }, label));
                });
                resultsSelect._records = data.results || [];
                resultsSelect.classList.toggle("hidden", (data.results || []).length === 0);
            } catch (err) { /* keep prior dropdown contents on transient search errors */ }
        }, 350);
    });

    resultsSelect.addEventListener("change", () => {
        if (resultsSelect.value === "") return;
        const rec = resultsSelect._records[parseInt(resultsSelect.value, 10)];
        dbCfg.taxonId = String(rec.context_id || "");
        dbCfg.taxonName = String(rec.context_name || "");
        refreshCaption();
        // Re-render just the run button's enabled state via a full panel pass
        // (cheap; this screen's data is small).
        render();
    });

    wrap.appendChild(input);
    wrap.appendChild(resultsSelect);
    wrap.appendChild(selectedCaption);
    return wrap;
}

function startOraAnalysisRun() {
    if (!_oraConfigIsRunnable()) return;
    const perDbApi = {};
    STATE.ora.analysisConfig.databases.forEach(db => {
        const cfg = STATE.ora.analysisConfig.perDatabase[db] || {};
        // Optional term-size filter (all databases). min defaults to 3; max blank -> no cap.
        const minTs = (cfg.minTermSize || cfg.minTermSize === 0) ? cfg.minTermSize : ORA_DEFAULT_MIN_TERM_SIZE;
        // "" is an explicit user opt-out (no cap); undefined (row never
        // rendered/touched) falls back to the same default the term-size
        // row itself uses, not to "no cap".
        const maxTs = cfg.maxTermSize === "" ? null : (typeof cfg.maxTermSize === "undefined" ? ORA_DEFAULT_MAX_TERM_SIZE : cfg.maxTermSize);
        if (db === ORA_CUSTOM_DB_LABEL) {
            perDbApi[db] = { id_type: cfg.idType || "accession", min_term_size: minTs, max_term_size: maxTs };
        } else {
            perDbApi[db] = {
                taxon_id: cfg.taxonId,
                taxon_name: cfg.taxonName,
                include_iea: typeof cfg.includeIea === "undefined" ? true : cfg.includeIea,
                min_term_size: minTs,
                max_term_size: maxTs,
            };
            // GO/Reactome only: flag whether the background checkbox is on AND a
            // file was actually uploaded, so the server keeps it (and prunes it
            // otherwise on this run) - unchecking disables a background even if
            // one is still stored locally from a previous upload.
            if (db === "GO" || db === "Reactome") {
                perDbApi[db].has_background = !!(cfg.useBackground && cfg.backgroundSummary);
            }
        }
    });
    const payload = {
        databases: STATE.ora.analysisConfig.databases,
        per_database: perDbApi,
    };
    apiPostJson("/api/ora/analysis/configure", payload).then(() => {
        STATE.ora.analysisRunLocked = true;
        STATE.ora.analysisResults = null;
        STATE.aaPage = "genomics_ora_results";
        render();
        return apiPostJson("/api/ora/analysis/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("ora-analysis-progress-inner");
                const label = document.getElementById("ora-analysis-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => {
                STATE.ora.analysisResults = result;
                render();
            },
            onError: (err) => {
                alert("ORA run failed: " + err);
                STATE.ora.analysisRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.ora.analysisRunLocked = false;
        render();
    });
}

function renderOraAnalysisResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Over-Representation Analysis: results"));

    if (!STATE.ora.analysisResults) {
        wrap.appendChild(alertBox("info", "Running analysis..."));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "ora-analysis-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "ora-analysis-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    Object.keys(STATE.ora.analysisResults).forEach((dbName) => {
        const dbResult = STATE.ora.analysisResults[dbName];
        if (dbName === "GO") {
            ["biological_process", "molecular_function", "cellular_component"].forEach((ns) => {
                if (dbResult[ns]) wrap.appendChild(renderOraResultBlock(dbName, ns, dbResult[ns]));
            });
        } else {
            const key = Object.keys(dbResult)[0];
            if (key) wrap.appendChild(renderOraResultBlock(dbName, null, dbResult[key]));
        }
    });

    const newRunBtn = el("button", { class: "btn" }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.ora.analysisConfig = { databases: [], perDatabase: {} };
        STATE.ora.analysisResults = null;
        STATE.ora.analysisRunLocked = false;
        STATE.aaPage = "genomics_ora_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

// ---------------------------------------------------------------------------
// Plotly chart helpers for ORA results
// ---------------------------------------------------------------------------

function _tableToObjects(table) {
    // Rows from df_payload(orient="records") are already dicts keyed by column name.
    return table.rows || [];
}

// sortBy: "fdr" (default, ascending - most significant first), "fold" (Fold
// Enrichment descending), "k" (gene count descending), "ratio" (k/n descending).
const ORA_SORT_COMPARATORS = {
    fdr: (a, b) => Number(a.FDR) - Number(b.FDR),
    fold: (a, b) => Number(b["Fold Enrichment"]) - Number(a["Fold Enrichment"]),
    k: (a, b) => Number(b.k) - Number(a.k),
    ratio: (a, b) => (Number(b.k) / Math.max(Number(b.n), 1)) - (Number(a.k) / Math.max(Number(a.n), 1)),
};

// What the chosen "Rank & show by" metric encodes on the dot-plot x-axis and the
// bar-plot bar length, so the ranked metric is actually visible (not just used to
// sort). For "fdr", the dot x stays Gene ratio because color already encodes FDR
// (avoids redundancy); the bar length is -log10(FDR) as before.
const ORA_METRIC_SPECS = {
    fdr: {
        barValue: r => -Math.log10(Number(r.FDR) + 1e-300), barTitle: "-log10(FDR)", barFmt: ".2f",
        dotValue: r => Number(r.k) / Math.max(Number(r.n), 1), dotTitle: "Gene ratio (k/n)", dotFmt: ".3f",
    },
    fold: {
        barValue: r => Number(r["Fold Enrichment"]), barTitle: "Fold Enrichment", barFmt: ".2f",
        dotValue: r => Number(r["Fold Enrichment"]), dotTitle: "Fold Enrichment", dotFmt: ".2f",
    },
    k: {
        barValue: r => Number(r.k), barTitle: "Gene count (k)", barFmt: ".0f",
        dotValue: r => Number(r.k), dotTitle: "Gene count (k)", dotFmt: ".0f",
    },
    ratio: {
        barValue: r => Number(r.k) / Math.max(Number(r.n), 1), barTitle: "Gene ratio (k/n)", barFmt: ".3f",
        dotValue: r => Number(r.k) / Math.max(Number(r.n), 1), dotTitle: "Gene ratio (k/n)", dotFmt: ".3f",
    },
};

// Shared, metric-agnostic hover block for dot/bar so the tooltip is always correct
// no matter what the x-axis encodes. customdata order must match: [name, fdr, fold, k, ratio].
const ORA_DOTBAR_HOVER = "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Fold enrichment: %{customdata[2]:.2f}<br>Gene count: %{customdata[3]}<br>Gene ratio: %{customdata[4]:.3f}<extra></extra>";

// alpha: significance threshold - rows with FDR >= alpha are excluded before
// sorting/slicing, so "Terms shown" can never pad a chart with non-significant
// terms. n === Infinity ("All significant") skips the cap entirely.
function _oraTopRows(table, n, sortBy, alpha) {
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : ORA_DEFAULT_SIG_THRESHOLD;
    const rows = _tableToObjects(table).filter(r =>
        r.FDR !== null && r.FDR !== "" && isFinite(Number(r.FDR)) && Number(r.FDR) < sigAlpha
    );
    const comparator = ORA_SORT_COMPARATORS[sortBy] || ORA_SORT_COMPARATORS.fdr;
    // Rows already arrive FDR-ascending from the server; re-sort only when a
    // different order is requested (sort is stable, so "fdr" is a no-op reorder).
    const sorted = sortBy && sortBy !== "fdr" ? rows.slice().sort(comparator) : rows;
    return n === Infinity ? sorted : sorted.slice(0, n);
}

function _oraChartId(kind, dbName, namespace) {
    return "chart-" + kind + "-" + dbName.replace(/\s/g, "_") + "-" + (namespace || "x");
}

function _escapePlotlyText(value) {
    return String(value || "")
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;");
}

function _wrapOraLabel(value, maxChars) {
    const words = String(value || "").split(/\s+/).filter(Boolean);
    const lines = [];
    let current = "";
    words.forEach((word) => {
        if (word.length > maxChars) {
            if (current) {
                lines.push(current);
                current = "";
            }
            for (let i = 0; i < word.length; i += maxChars) {
                lines.push(word.slice(i, i + maxChars));
            }
            return;
        }
        if (!current) {
            current = word;
        } else if ((current + " " + word).length <= maxChars) {
            current += " " + word;
        } else {
            lines.push(current);
            current = word;
        }
    });
    if (current) lines.push(current);
    return (lines.length ? lines : [""]).map(_escapePlotlyText).join("<br>");
}

// Height grows (Y axis only) with the number of terms AND the chosen font size,
// so labels never overlap regardless of how many terms are shown or how large the
// text is. Plotly spaces categories evenly, so each term's slot must be tall enough
// for the TALLEST wrapped label at that font size; total height = terms x slot +
// margins (min floor for small sets). Must be recomputed on every rebuild (top-N /
// sort / font-size change) - it is never cached - so the guarantee always holds.
function _oraCategoricalChartHeight(wrappedLabels, fontSize) {
    const size = Number(fontSize) > 0 ? Number(fontSize) : 12;
    const LINE_PX = Math.max(16, Math.round(size * 1.7));   // vertical space per line of label text, scaled to font size
    const ROW_MIN = Math.max(28, LINE_PX + 8);               // minimum comfortable height per term slot
    const MARGINS = 22 + 56;    // top + bottom margins (must match the layout below)
    const MIN_HEIGHT = 360;     // floor so a handful of terms isn't a giant chart
    const n = Math.max(1, wrappedLabels.length);
    const maxLines = wrappedLabels.reduce((mx, label) => {
        return Math.max(mx, String(label).split("<br>").length);
    }, 1);
    const perSlot = Math.max(ROW_MIN, maxLines * LINE_PX);
    return Math.max(MIN_HEIGHT, n * perSlot + MARGINS);
}

// Marker diameters for a count-encoded dot/bubble chart.
//
// Plotly's marker.size is a raw pixel DIAMETER with nothing applied to it
// (no sizemode/sizeref anywhere in this file), so the old per-chart formulas
// - `Math.max(7, k * 4)` and friends - scaled without bound: an ORA term
// holding k=200 of the user's genes produced an 800px circle on a ~900px
// plot, i.e. one point covering the chart. Every test fixture was small
// enough (k <= 19) to hide it, and no test asserted on marker size at all.
//
// Two fixes in one:
//   1. the largest value always maps to maxPx REGARDLESS of dataset scale,
//      so a 19-gene fixture and a 20,000-gene list both render sensibly;
//   2. sqrt() puts the value on AREA rather than diameter, which is how
//      magnitude is actually perceived - the old diameter-proportional
//      encoding overstated big terms quadratically (4x the genes drew 16x
//      larger), so this is a correctness fix, not only a cosmetic one.
function _magnitudeMarkerSizes(values, minPx, maxPx) {
    const lo = Number(minPx) > 0 ? Number(minPx) : 8;
    const hi = Number(maxPx) > 0 ? Number(maxPx) : 30;
    const vals = values.map((v) => {
        const n = Number(v);
        return isFinite(n) && n > 0 ? n : 0;
    });
    // reduce(), NOT Math.max(...vals): "All significant" can pass thousands
    // of terms, and spreading an array that large risks a stack overflow.
    let maxV = 0;
    let minV = Infinity;
    for (let i = 0; i < vals.length; i++) {
        if (vals[i] > maxV) maxV = vals[i];
        if (vals[i] < minV) minV = vals[i];
    }
    if (!(maxV > 0)) return vals.map(() => lo);
    // Every term the same count: a chart of uniformly maximal circles would
    // read as "everything is huge" when it actually means "no variation".
    if (minV === maxV) return vals.map(() => (lo + hi) / 2);
    return vals.map((v) => lo + Math.sqrt(v / maxV) * (hi - lo));
}

function _oraCameraIcon() {
    if (typeof Plotly !== "undefined" && Plotly.Icons && Plotly.Icons.camera) return Plotly.Icons.camera;
    return {
        width: 1000,
        ascent: 1000,
        descent: 0,
        path: "M120 260h180l55-80h290l55 80h180c45 0 80 35 80 80v400c0 45-35 80-80 80H120c-45 0-80-35-80-80V340c0-45 35-80 80-80zm380 460c105 0 190-85 190-190s-85-190-190-190-190 85-190 190 85 190 190 190zm0-70c-66 0-120-54-120-120s54-120 120-120 120 54 120 120-54 120-120 120z",
    };
}

async function _downloadOraChart(div, dbName, namespace, kindLabel, kind, fmt) {
    if (typeof Plotly === "undefined") { alert("Plotly not loaded."); return; }
    if (!div || !div.data) {
        alert(kindLabel + " chart is not rendered yet. Show it first, then export.");
        return;
    }
    const filename = dbName + (namespace ? "_" + namespace : "") + "_" + kind;
    const fullLayout = div._fullLayout || {};
    const height = Math.max(500, Math.ceil(fullLayout.height || div.clientHeight || 500));

    if (fmt === "tiff") {
        // Plotly.js cannot export TIFF directly - render a high-res PNG client-side,
        // then ask the server to re-encode it as a 300 DPI TIFF via Pillow.
        try {
            const scale = 300 / 96; // scale a 96-dpi-equivalent render up to 300 dpi
            const dataUrl = await Plotly.toImage(div, { format: "png", width: 900, height, scale });
            const resp = await apiFetch("/api/ora/chart-tiff", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ image: dataUrl, filename: filename + ".tiff" }),
            });
            const blob = await resp.blob();
            const url = URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url; a.download = filename + ".tiff";
            document.body.appendChild(a); a.click(); a.remove();
            URL.revokeObjectURL(url);
        } catch (err) {
            alert("TIFF export failed: " + (err.message || err));
        }
        return;
    }

    Plotly.downloadImage(div, { format: fmt, filename: filename, width: 900, height });
}

function _showOraChartExportMenu(chartDiv, dbName, namespace, kind, kindLabel) {
    const existing = document.querySelector(".ora-chart-export-menu");
    if (existing) existing.remove();

    const menu = el("div", { class: "ora-chart-export-menu" });
    const fmts = [["SVG", "svg"], ["PNG", "png"], ["TIFF (300 dpi)", "tiff"]];
    fmts.forEach(([fmtLabel, fmt]) => {
        const item = el("button", { type: "button", class: "ora-chart-export-item" }, kindLabel + " " + fmtLabel);
        item.addEventListener("click", () => {
            menu.remove();
            _downloadOraChart(chartDiv, dbName, namespace, kindLabel, kind, fmt);
        });
        menu.appendChild(item);
    });
    document.body.appendChild(menu);

    const rect = chartDiv.getBoundingClientRect();
    menu.style.top = Math.max(8, rect.top + window.scrollY + 28) + "px";
    menu.style.left = Math.max(8, rect.right + window.scrollX - menu.offsetWidth - 10) + "px";

    setTimeout(() => {
        document.addEventListener("click", function closeMenu(evt) {
            if (!menu.contains(evt.target)) {
                menu.remove();
                document.removeEventListener("click", closeMenu);
            }
        });
    }, 0);
}

function _oraChartConfig(dbName, namespace, kind, kindLabel) {
    return {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ["toImage"],
        modeBarButtonsToAdd: [{
            name: "Download " + (kindLabel || "chart"),
            title: "Download " + (kindLabel || "chart"),
            icon: _oraCameraIcon(),
            click: (chartDiv) => _showOraChartExportMenu(chartDiv, dbName, namespace, kind, kindLabel || "Chart"),
        }],
    };
}

function _clearOraChartSelection(div) {
    if (typeof Plotly === "undefined" || !div) return;
    try {
        Plotly.relayout(div, { dragmode: "zoom", selections: [] });
        if (div.data && div.data.length) {
            Plotly.restyle(div, { selectedpoints: null });
        }
    } catch (e) { /* selection state is best-effort UI cleanup */ }
}

function _enableOraSelectToolToggle(div) {
    if (!div || div._oraSelectToolToggleBound) return;
    div._oraSelectToolToggleBound = true;
    div.addEventListener("click", (evt) => {
        const button = evt.target.closest ? evt.target.closest(".modebar-btn") : null;
        if (!button) return;
        const title = (
            button.getAttribute("data-title") ||
            button.getAttribute("aria-label") ||
            button.getAttribute("title") ||
            ""
        ).toLowerCase();
        const isBoxSelect = title.indexOf("box select") !== -1 || title === "select";
        const isLassoSelect = title.indexOf("lasso select") !== -1 || title === "lasso";
        if (!isBoxSelect && !isLassoSelect) return;
        const dragmode = div._fullLayout && div._fullLayout.dragmode;
        const wasActiveSelectTool = (isBoxSelect && dragmode === "select") || (isLassoSelect && dragmode === "lasso");
        if (wasActiveSelectTool) {
            setTimeout(() => _clearOraChartSelection(div), 0);
        }
    }, true);
}

function _newOraPlot(div, traces, layout, dbName, namespace, kind, kindLabel, onReady) {
    const config = _oraChartConfig(dbName, namespace, kind, kindLabel);
    function draw() {
        const plot = Plotly.newPlot(div, traces, layout, config);
        if (plot && typeof plot.then === "function") {
            plot.then(() => {
                _enableOraSelectToolToggle(div);
                if (onReady) onReady();
            });
        } else {
            _enableOraSelectToolToggle(div);
            if (onReady) onReady();
        }
    }
    // Stashed so a chart first drawn while its result-nav pane was still
    // off-DOM/display:none (width 0) can be fully REBUILT once revealed.
    // Confirmed live: Plotly.Plots.resize() - the fix makeResultNav's
    // activate() used to rely on exclusively - rescales the existing
    // layout but does NOT recompute things like legend height/position
    // that were miscalculated against the wrong original width (measured:
    // the legend stayed pinned ~20px too high, overlapping the xaxis title,
    // even after resize()/relayout({})/Plotly.react() - only a genuinely
    // fresh newPlot() at the real container size fixes it). See activate().
    div._omiscerRedraw = draw;
    draw();
}

function buildDotPlot(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const rows = _oraTopRows(table, topN, opts.sortBy, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const spec = ORA_METRIC_SPECS[opts.sortBy] || ORA_METRIC_SPECS.fdr;
    const fdrVals = rows.map(r => Number(r.FDR));
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const trace = {
        // x = the chosen "Rank & show by" metric (falls back to gene ratio for FDR).
        x: rows.map(spec.dotValue),
        y: wrappedNames,
        // customdata = [name, FDR, fold enrichment, gene count, gene ratio] so the
        // hover is complete and correct regardless of what the x-axis encodes.
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.FDR), Number(r["Fold Enrichment"]),
            Number(r.k), Number(r.k) / Math.max(Number(r.n), 1),
        ]),
        mode: "markers",
        type: "scatter",
        marker: {
            size: _magnitudeMarkerSizes(rows.map(r => Number(r.k)), 8, 30),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.7 },
        },
        hovertemplate: ORA_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: spec.dotTitle, font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: false },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 105, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "dot", "Dot plot");
}

function buildBarPlot(div, table, dbName, namespace, barColor, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const rows = _oraTopRows(table, topN, opts.sortBy, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    const color = ORA_BAR_COLORS.some((item) => item[1] === barColor) ? barColor : ORA_BAR_COLORS[0][1];
    const spec = ORA_METRIC_SPECS[opts.sortBy] || ORA_METRIC_SPECS.fdr;
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const trace = {
        // Bar length = the chosen "Rank & show by" metric; rows are already sorted by
        // it (via _oraTopRows), so with y reversed the bars descend from the top.
        x: rows.map(spec.barValue),
        y: wrappedNames,
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.FDR), Number(r["Fold Enrichment"]),
            Number(r.k), Number(r.k) / Math.max(Number(r.n), 1),
        ]),
        type: "bar",
        orientation: "h",
        marker: { color },
        hovertemplate: ORA_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: spec.barTitle, font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: false },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 70, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "bar", "Bar plot");
}

function buildVolcanoPlot(div, table, dbName, namespace, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : 0.05;
    const rows = _tableToObjects(table).filter(r => r["Fold Enrichment"] !== null && isFinite(Number(r["Fold Enrichment"])));
    if (!rows.length) { div.textContent = "No terms to plot."; return; }
    const logFE = rows.map(r => Math.log2(Math.max(Number(r["Fold Enrichment"]), 0.001)));
    // Y axis: -log10(p-value) by default, or -log10(FDR) for consistency with the
    // FDR-based point coloring (opts.volcanoY === "fdr").
    const useFdrY = opts.volcanoY === "fdr";
    const yField = useFdrY ? "FDR" : "p-value";
    const yTitle = useFdrY ? "-log10(FDR)" : "-log10(p-value)";
    const yHoverLabel = useFdrY ? "-log10(FDR)" : "-log10(p)";
    const logY = rows.map(r => -Math.log10(Number(r[yField]) + 1e-300));
    const colors = rows.map(r => Number(r.FDR) < sigAlpha ? "#c0392b" : "#95a5a6");
    const trace = {
        x: logFE,
        y: logY,
        mode: "markers",
        type: "scatter",
        text: rows.map(r => r["Term Name"]),
        marker: { color: colors, size: 8, opacity: 0.8 },
        hovertemplate: "%{text}<br>log2(FE)=%{x:.2f}<br>" + yHoverLabel + "=%{y:.2f}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "log2(Fold Enrichment)", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { title: { text: yTitle, font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: false },
        height: 500,
        margin: { l: 95, r: 110, t: 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "volcano", "Volcano");
}

// ---------------------------------------------------------------------------
// Gene-set Burden Test charts (dot/bar/volcano). Burden's result columns
// differ from ORA's (Genes in set / Effect (mean difference) / Mean score,
// rather than k/K/n/N/Fold Enrichment), so these are separate builders - but
// they reuse every generic piece above (_newOraPlot, _oraChartId,
// _oraCategoricalChartHeight, _wrapOraLabel) so charts look and export
// identically to ORA's.
// ---------------------------------------------------------------------------

// Burden's dot plot has a single, fixed encoding (no "Rank & show by" toggle -
// unlike ORA, burden only has one other continuous metric besides FDR, so a
// metric-switcher control had nothing meaningful to switch between on the dot
// plot: x-axis is always Effect (mean difference)). The bar plot DOES offer a
// metric choice (see "Bar metric" in renderBurdenResultBlock) - it's scoped to
// just the bar plot, so it doesn't repeat the earlier dot-plot confusion.
const BURDEN_DOTBAR_HOVER = "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Effect (mean diff): %{customdata[2]:.3f}<br>Genes in set: %{customdata[3]}<extra></extra>";

// alpha: significance threshold - rows with FDR >= alpha are excluded before
// slicing (same semantics as _oraTopRows). n === Infinity ("All significant")
// skips the cap entirely. sortBy: "fdr" (default; the table already arrives
// FDR-ascending from the server, so this is a no-op reorder) or "effect"
// (re-sorts by Effect (mean difference) descending, e.g. for the bar plot's
// "Bar metric" = Effect mode, so its bars still come out as a clean ranked
// staircase rather than jumping around in FDR order).
function _burdenTopRows(table, n, alpha, sortBy) {
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : ORA_DEFAULT_SIG_THRESHOLD;
    const rows = _tableToObjects(table).filter(r =>
        r.FDR !== null && r.FDR !== "" && isFinite(Number(r.FDR)) && Number(r.FDR) < sigAlpha
    );
    const sorted = sortBy === "effect"
        ? rows.slice().sort((a, b) => Number(b["Effect (mean difference)"]) - Number(a["Effect (mean difference)"]))
        : rows;
    return n === Infinity ? sorted : sorted.slice(0, n);
}

function buildBurdenDotPlot(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const rows = _burdenTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const fdrVals = rows.map(r => Number(r.FDR));
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const trace = {
        x: rows.map(r => Number(r["Effect (mean difference)"])),
        y: wrappedNames,
        // customdata = [name, FDR, effect, genes-in-set] so hover is complete
        // regardless of what the x-axis encodes.
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.FDR), Number(r["Effect (mean difference)"]),
            Number(r["Genes in set"]),
        ]),
        mode: "markers",
        type: "scatter",
        marker: {
            size: _magnitudeMarkerSizes(rows.map(r => Number(r["Genes in set"])), 8, 30),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.7 },
        },
        hovertemplate: BURDEN_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Effect (mean difference)", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 105, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "dot", "Dot plot");
}

function buildBurdenBarPlot(div, table, dbName, namespace, barColor, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    // "Bar metric" ("fdr" default | "effect") picks both the bar length AND the
    // sort order, so the bars always come out as a clean descending staircase
    // regardless of which metric is chosen (see _burdenTopRows).
    const barMetric = opts.barMetric === "effect" ? "effect" : "fdr";
    const rows = _burdenTopRows(table, topN, opts.alpha, barMetric);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    const color = ORA_BAR_COLORS.some((item) => item[1] === barColor) ? barColor : ORA_BAR_COLORS[0][1];
    const names = rows.map(r => String(r["Term Name"] || ""));
    const wrappedNames = names.map(name => _wrapOraLabel(name, 44));
    const chartHeight = _oraCategoricalChartHeight(wrappedNames, fontSize);
    const barValue = barMetric === "effect"
        ? (r => Number(r["Effect (mean difference)"]))
        : (r => -Math.log10(Number(r.FDR) + 1e-300));
    const barTitle = barMetric === "effect" ? "Effect (mean difference)" : "-log10(FDR)";
    const trace = {
        x: rows.map(barValue),
        y: wrappedNames,
        customdata: rows.map(r => [
            String(r["Term Name"] || ""), Number(r.FDR), Number(r["Effect (mean difference)"]),
            Number(r["Genes in set"]),
        ]),
        type: "bar",
        orientation: "h",
        marker: { color },
        hovertemplate: BURDEN_DOTBAR_HOVER,
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: barTitle, font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { autorange: "reversed", tickfont: { size: fontSize }, automargin: true },
        height: chartHeight,
        margin: { l: 385, r: 70, t: 22, b: 56 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "bar", "Bar plot");
}

function buildBurdenVolcanoPlot(div, table, dbName, namespace, alpha, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const sigAlpha = Number(alpha) > 0 ? Number(alpha) : 0.05;
    const rows = _tableToObjects(table).filter(r => r["Effect (mean difference)"] !== null && isFinite(Number(r["Effect (mean difference)"])));
    if (!rows.length) { div.textContent = "No terms to plot."; return; }
    // x = Effect (mean difference) directly - unlike ORA's fold enrichment (a
    // ratio, so plotted log2), burden's effect is already a linear difference
    // of means, so no log transform is needed.
    const effect = rows.map(r => Number(r["Effect (mean difference)"]));
    const useFdrY = opts.volcanoY === "fdr";
    const yField = useFdrY ? "FDR" : "p-value";
    const yTitle = useFdrY ? "-log10(FDR)" : "-log10(p-value)";
    const yHoverLabel = useFdrY ? "-log10(FDR)" : "-log10(p)";
    const logY = rows.map(r => -Math.log10(Number(r[yField]) + 1e-300));
    const colors = rows.map(r => Number(r.FDR) < sigAlpha ? "#c0392b" : "#95a5a6");
    const trace = {
        x: effect,
        y: logY,
        mode: "markers",
        type: "scatter",
        text: rows.map(r => r["Term Name"]),
        marker: { color: colors, size: 8, opacity: 0.8 },
        hovertemplate: "%{text}<br>Effect=%{x:.3f}<br>" + yHoverLabel + "=%{y:.2f}<extra></extra>",
    };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Effect (mean difference)", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { title: { text: yTitle, font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: false },
        height: 500,
        margin: { l: 95, r: 110, t: 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, dbName, namespace, "volcano", "Volcano");
}

// ---------------------------------------------------------------------------
// Relational charts: enrichment map (terms linked by shared genes) and
// gene-concept network (terms <-> their overlapping genes). Both are laid out
// with a small seeded force-directed algorithm and rendered as plain Plotly
// scatter traces, so they get the same interactive toolbar + SVG/PNG/TIFF
// export as every other ORA chart with zero new dependencies.
// ---------------------------------------------------------------------------

function _oraGeneSet(row) {
    // Proteomics names this column "Overlapping Proteins" and formats each
    // member as "ACCESSION (SYMBOL)" (see Proteomics_module.OVERLAP_COLUMN).
    // Read either column name, and keep only the leading identifier - the
    // callers that consume this set (gene/protein-concept network node ids,
    // and its log2FC lookup) key on the bare accession/symbol, so passing
    // "P05386 (RPLP0)" through would break every lookup silently. Genomics
    // and Transcriptomics values contain no parenthesis, so the strip is a
    // no-op for them.
    const raw = String(row["Overlapping Genes"] || row["Overlapping Proteins"] || "").trim();
    if (!raw) return new Set();
    return new Set(
        raw.split(",")
           .map(g => g.trim().split(" (")[0].trim())
           .filter(Boolean)
    );
}

// Burden analog of _oraGeneSet: the burden result table's per-term gene list
// lives in the "Genes" column (the scored genes in that term), not "Overlapping
// Genes".
function _burdenGeneSet(row) {
    const raw = String(row["Genes"] || "").trim();
    if (!raw) return new Set();
    return new Set(raw.split(",").map(g => g.trim()).filter(Boolean));
}

// Minimal seeded Fruchterman-Reingold force layout, adequate for the small
// graphs here (<=~20 terms, or terms + their handful of overlapping genes).
// nodeIds: array of unique string ids. edges: array of {source, target, weight}.
// Returns {id: {x, y}} with coordinates roughly in [-1, 1] x [-1, 1].
// Seeded (mulberry32) so the same graph lays out identically across re-renders
// (e.g. toggling a color control shouldn't scramble node positions).
function _forceLayout(nodeIds, edges, opts) {
    opts = opts || {};
    const n = nodeIds.length;
    const pos = {};
    if (n === 0) return pos;

    const iterations = opts.iterations || 150;
    const area = 4;
    const k = Math.sqrt(area / n);

    let seed = opts.seed || 42;
    function rand() {
        seed = (seed + 0x6D2B79F5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }

    nodeIds.forEach((id) => {
        pos[id] = { x: rand() * 2 - 1, y: rand() * 2 - 1 };
    });
    if (n === 1) return pos;

    const idxOf = {};
    nodeIds.forEach((id, i) => { idxOf[id] = i; });
    const validEdges = (edges || []).filter(e => idxOf[e.source] !== undefined && idxOf[e.target] !== undefined);

    let temperature = 0.1;
    const cooling = temperature / iterations;

    for (let iter = 0; iter < iterations; iter++) {
        const disp = nodeIds.map(() => ({ x: 0, y: 0 }));

        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                let dx = pos[nodeIds[i]].x - pos[nodeIds[j]].x;
                let dy = pos[nodeIds[i]].y - pos[nodeIds[j]].y;
                const dist = Math.sqrt(dx * dx + dy * dy) || 0.001;
                const force = (k * k) / dist;
                dx = (dx / dist) * force; dy = (dy / dist) * force;
                disp[i].x += dx; disp[i].y += dy;
                disp[j].x -= dx; disp[j].y -= dy;
            }
        }

        validEdges.forEach((e) => {
            const i = idxOf[e.source], j = idxOf[e.target];
            let dx = pos[nodeIds[i]].x - pos[nodeIds[j]].x;
            let dy = pos[nodeIds[i]].y - pos[nodeIds[j]].y;
            const dist = Math.sqrt(dx * dx + dy * dy) || 0.001;
            const force = ((dist * dist) / k) * (e.weight || 1);
            dx = (dx / dist) * force; dy = (dy / dist) * force;
            disp[i].x -= dx; disp[i].y -= dy;
            disp[j].x += dx; disp[j].y += dy;
        });

        for (let i = 0; i < n; i++) {
            const dLen = Math.sqrt(disp[i].x * disp[i].x + disp[i].y * disp[i].y) || 0.001;
            const capped = Math.min(dLen, temperature);
            pos[nodeIds[i]].x += (disp[i].x / dLen) * capped;
            pos[nodeIds[i]].y += (disp[i].y / dLen) * capped;
        }
        temperature = Math.max(0.001, temperature - cooling);
    }

    return pos;
}

// Shared network renderer: hides axes/grid (positions are layout coordinates,
// not data), keeps a 1:1 aspect ratio, and routes through _newOraPlot so
// network charts get the standard toolbar + per-chart export menu.
function _plotNetwork(div, traces, layoutOverrides, dbName, namespace, kind, kindLabel, onReady) {
    if (typeof Plotly === "undefined") return;
    const layout = Object.assign({
        showlegend: false,
        hovermode: "closest",
        // Default click-drag is disabled (rather than Plotly's box-zoom) so a plain
        // click-drag moves a node instead; zoom/pan remain available via the modebar.
        dragmode: false,
        xaxis: { visible: false, zeroline: false, showgrid: false },
        yaxis: { visible: false, zeroline: false, showgrid: false, scaleanchor: "x", scaleratio: 1 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
        margin: { l: 20, r: 20, t: 20, b: 20 },
        height: 620,
    }, layoutOverrides || {});
    _newOraPlot(div, traces, layout, dbName, namespace, kind, kindLabel, onReady);
}

// Makes term/gene nodes on a network chart individually draggable by click-drag.
// Plotly has no built-in "draggable scatter point" interaction, so this implements
// one: hit-test against each node's current PIXEL position (derived each call from
// the chart's own axis ranges + plot-area geometry, both stable once rendered),
// then on drag directly mutate the node/edge traces' x/y arrays and Plotly.redraw.
// Requires dragmode:false (set in _plotNetwork) so box-zoom doesn't fight with it;
// modebar Zoom/Pan/Reset buttons are unaffected and still work.
//   div          - the rendered chart div (must already have _fullLayout, i.e.
//                  called from _newOraPlot's onReady, after Plotly.newPlot resolves)
//   traceGroups  - [{traceIndex, nodeIds}], one entry per draggable trace; nodeIds
//                  must be in the same order as that trace's x/y arrays
//   pos          - the live {nodeId: {x,y}} map from _forceLayout (mutated in place)
//   edgeTraceIndex - index of the single edges line-trace sharing the same pos space
//   edges        - the edge list used to build the edge trace ({source,target,...})
function _makeNetworkDraggable(div, traceGroups, pos, edgeTraceIndex, edges) {
    if (typeof Plotly === "undefined" || !div) return;

    const nodeLocation = {};
    traceGroups.forEach(({ traceIndex, nodeIds }) => {
        nodeIds.forEach((id, i) => { nodeLocation[id] = { traceIndex, localIndex: i }; });
    });

    const edgeSlotsByNode = {};
    (edges || []).forEach((e, idx) => {
        (edgeSlotsByNode[e.source] = edgeSlotsByNode[e.source] || []).push({ idx, role: "source" });
        (edgeSlotsByNode[e.target] = edgeSlotsByNode[e.target] || []).push({ idx, role: "target" });
    });

    function pixelToData(clientX, clientY) {
        const fl = div._fullLayout;
        const rect = div.getBoundingClientRect();
        const size = fl._size;
        const xr = fl.xaxis.range, yr = fl.yaxis.range;
        const px = clientX - rect.left, py = clientY - rect.top;
        return {
            x: xr[0] + ((px - size.l) / size.w) * (xr[1] - xr[0]),
            y: yr[1] - ((py - size.t) / size.h) * (yr[1] - yr[0]),
        };
    }

    function dataToPixel(x, y) {
        const fl = div._fullLayout;
        const size = fl._size;
        const xr = fl.xaxis.range, yr = fl.yaxis.range;
        return {
            px: size.l + ((x - xr[0]) / (xr[1] - xr[0])) * size.w,
            py: size.t + ((yr[1] - y) / (yr[1] - yr[0])) * size.h,
        };
    }

    function hitTest(clientX, clientY) {
        const rect = div.getBoundingClientRect();
        const mx = clientX - rect.left, my = clientY - rect.top;
        const TOL = 20; // pixels
        let best = null, bestDist = TOL;
        Object.keys(nodeLocation).forEach((id) => {
            const p = pos[id];
            if (!p) return;
            const { px, py } = dataToPixel(p.x, p.y);
            const d = Math.sqrt((px - mx) * (px - mx) + (py - my) * (py - my));
            if (d < bestDist) { bestDist = d; best = id; }
        });
        return best;
    }

    function applyMove(id) {
        const loc = nodeLocation[id];
        if (loc && div.data[loc.traceIndex]) {
            div.data[loc.traceIndex].x[loc.localIndex] = pos[id].x;
            div.data[loc.traceIndex].y[loc.localIndex] = pos[id].y;
        }
        const slots = edgeSlotsByNode[id] || [];
        const edgeTrace = div.data[edgeTraceIndex];
        if (edgeTrace) {
            slots.forEach(({ idx, role }) => {
                const base = idx * 3 + (role === "source" ? 0 : 1);
                edgeTrace.x[base] = pos[id].x;
                edgeTrace.y[base] = pos[id].y;
            });
        }
    }

    let draggingId = null;
    let rafPending = false;
    function scheduleRedraw() {
        if (rafPending) return;
        rafPending = true;
        requestAnimationFrame(() => {
            rafPending = false;
            if (document.body.contains(div)) Plotly.redraw(div);
        });
    }

    function onMouseDown(evt) {
        if (evt.button !== 0) return;
        const hit = hitTest(evt.clientX, evt.clientY);
        if (!hit) return;
        draggingId = hit;
        evt.preventDefault();
        evt.stopPropagation();
    }
    function onMouseMove(evt) {
        if (!draggingId) return;
        pos[draggingId] = pixelToData(evt.clientX, evt.clientY);
        applyMove(draggingId);
        scheduleRedraw();
        evt.preventDefault();
    }
    function onMouseUp() { draggingId = null; }

    div.addEventListener("mousedown", onMouseDown, true);
    document.addEventListener("mousemove", onMouseMove);
    document.addEventListener("mouseup", onMouseUp);
}

function buildEnrichmentMap(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const edgeThreshold = opts.edgeThreshold || 1;
    const rows = _oraTopRows(table, topN, opts.sortBy, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    if (rows.length < 2) { div.textContent = "Need at least 2 terms to draw an enrichment map."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const geneSets = rows.map(_oraGeneSet);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r.FDR));
    const kVals = rows.map(r => Number(r.k));
    const foldVals = rows.map(r => Number(r["Fold Enrichment"]));

    const edges = [];
    for (let i = 0; i < rows.length; i++) {
        for (let j = i + 1; j < rows.length; j++) {
            let overlap = 0;
            geneSets[i].forEach((g) => { if (geneSets[j].has(g)) overlap++; });
            if (overlap < edgeThreshold) continue;
            const union = geneSets[i].size + geneSets[j].size - overlap;
            const jaccard = union > 0 ? overlap / union : 0;
            edges.push({ source: termIds[i], target: termIds[j], weight: jaccard, overlap });
        }
    }

    const pos = _forceLayout(termIds, edges, { iterations: 150 });

    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#B7C6D6", width: 1.5 }, hoverinfo: "skip", showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const nodeTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], fdrVals[i], kVals[i], foldVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(kVals, 10, 34),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.7 },
            line: { color: "#2c3e50", width: 1 },
        },
        hovertemplate: "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Genes: %{customdata[2]}<br>Fold enrichment: %{customdata[3]:.2f}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, nodeTrace], {}, dbName, namespace, "enrichmap", "Enrichment map",
        () => _makeNetworkDraggable(div, [{ traceIndex: 1, nodeIds: termIds }], pos, 0, edges));
}

function buildGeneConceptNetwork(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    // Gene-concept networks get crowded fast (every gene becomes its own node),
    // so cap top-N tighter than the dot/bar default even if a larger value is set.
    const topN = opts.topN === Infinity ? 15 : Math.min(opts.topN || 20, 15);
    const rows = _oraTopRows(table, topN, opts.sortBy, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r.FDR));
    const kVals = rows.map(r => Number(r.k));
    const foldVals = rows.map(r => Number(r["Fold Enrichment"]));
    const geneSets = rows.map(_oraGeneSet);

    const geneIndex = {};
    const geneIds = [];
    geneSets.forEach((set) => {
        set.forEach((g) => {
            if (!(g in geneIndex)) { geneIndex[g] = "G" + geneIds.length; geneIds.push(g); }
        });
    });

    const nodeIds = termIds.concat(geneIds.map(g => geneIndex[g]));
    const edges = [];
    rows.forEach((r, i) => {
        geneSets[i].forEach((g) => edges.push({ source: termIds[i], target: geneIndex[g], weight: 1 }));
    });

    const pos = _forceLayout(nodeIds, edges, { iterations: 150 });

    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#D8E1EA", width: 1 }, hoverinfo: "skip", showlegend: false,
    };

    const geneTrace = {
        x: geneIds.map(g => pos[geneIndex[g]].x),
        y: geneIds.map(g => pos[geneIndex[g]].y),
        mode: "markers",
        type: "scatter",
        text: geneIds,
        hovertemplate: "%{text}<extra></extra>",
        marker: { size: 8, color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } },
        showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const termTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], fdrVals[i], kVals[i], foldVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(kVals, 14, 38),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.7 },
            line: { color: "#2c3e50", width: 1.5 },
            symbol: "diamond",
        },
        hovertemplate: "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Genes: %{customdata[2]}<br>Fold enrichment: %{customdata[3]:.2f}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, geneTrace, termTrace], {}, dbName, namespace, "network", "Gene-concept network",
        () => _makeNetworkDraggable(
            div,
            [{ traceIndex: 1, nodeIds: geneIds.map(g => geneIndex[g]) }, { traceIndex: 2, nodeIds: termIds }],
            pos, 0, edges
        ));
}

// ---------------------------------------------------------------------------
// Transcriptomics ORA — two transcriptomics-specific charts. Both reuse the
// generic plumbing above (_newOraPlot/_plotNetwork, _forceLayout, _wrapOraLabel,
// _oraGeneSet, ORA_DOT_COLOR_SCALES) but read per-GENE log2FC/padj, which
// Genomics' ORA has no equivalent of (it has no fold-change axis).
// ---------------------------------------------------------------------------

// Gene-level volcano: one point per TESTED gene (not per term) - the standard
// RNA-seq figure, log2FC vs -log10(padj), colored by up/down/not-significant
// against this run's own thresholds. Replaces Genomics' term-level volcano
// (fold-enrichment vs FDR per GO term), which is unusual for transcriptomics
// and is intentionally not offered here.
function buildTranscriptomicsGeneVolcano(div, resolvedRows, slug, padjCutoff, lfcCutoff, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const rows = (resolvedRows || []).filter(r =>
        r.log2FC !== null && isFinite(Number(r.log2FC)) && r.padj !== null && isFinite(Number(r.padj)) && Number(r.padj) > 0
    );
    if (!rows.length) { div.textContent = "No genes with a usable log2FC/padj to plot."; return; }

    const x = rows.map(r => Number(r.log2FC));
    const y = rows.map(r => -Math.log10(Number(r.padj) + 1e-300));
    const labels = rows.map(r => String(r.gene_symbol || r.input_id || ""));
    const status = rows.map((r) => {
        const isSig = Number(r.padj) <= padjCutoff && Math.abs(Number(r.log2FC)) >= lfcCutoff;
        if (!isSig) return "ns";
        return Number(r.log2FC) >= lfcCutoff ? "up" : "down";
    });
    const colorByStatus = { up: "#B84A52", down: "#1B639C", ns: "#B7C0C9" };
    const colors = status.map(s => colorByStatus[s]);

    const trace = {
        x, y, mode: "markers", type: "scatter",
        text: labels,
        customdata: rows.map((r, i) => [labels[i], Number(r.log2FC), Number(r.padj), status[i]]),
        marker: { color: colors, size: 7, opacity: 0.75, line: { color: "#2c3e50", width: 0.3 } },
        hovertemplate: "%{customdata[0]}<br>log2FC=%{customdata[1]:.2f}<br>padj=%{customdata[2]:.3e}<br>%{customdata[3]}<extra></extra>",
    };

    const maxAbsX = Math.max(1, ...x.map(v => Math.abs(v)));
    const maxY = Math.max(1, ...y);
    const shapes = [
        { type: "line", x0: lfcCutoff, x1: lfcCutoff, y0: 0, y1: maxY * 1.05, line: { color: "#999", width: 1, dash: "dash" } },
        { type: "line", x0: -lfcCutoff, x1: -lfcCutoff, y0: 0, y1: maxY * 1.05, line: { color: "#999", width: 1, dash: "dash" } },
        { type: "line", x0: -maxAbsX * 1.1, x1: maxAbsX * 1.1, y0: -Math.log10(padjCutoff), y1: -Math.log10(padjCutoff), line: { color: "#999", width: 1, dash: "dash" } },
    ];

    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "log2 Fold-Change", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: true, zerolinecolor: "#ccc" },
        yaxis: { title: { text: "-log10(padj)", font: { size: fontSize + 4 } }, tickfont: { size: fontSize + 1 }, zeroline: false },
        shapes,
        height: 500,
        margin: { l: 95, r: 40, t: 24, b: 64 },
        plot_bgcolor: "#FAFCFE",
        paper_bgcolor: "#FAFCFE",
    };
    _newOraPlot(div, [trace], layout, slug, null, "genevolcano", "Gene volcano");
}

// Gene-concept network colored by log2FC (standard cnetplot(foldChange=...)
// practice): term/category nodes keep the FDR-significance color+size encoding
// (as in buildGeneConceptNetwork), but gene nodes are colored on a diverging
// red-blue scale centered at 0 by their log2FC, so it's visible at a glance
// which genes driving a term went up vs down.
function buildTranscriptomicsGeneConceptNetwork(div, table, slug, geneLog2FC, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const topN = opts.topN === Infinity ? 15 : Math.min(opts.topN || 20, 15);
    const rows = _oraTopRows(table, topN, opts.sortBy, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === opts.dotColorScale) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r.FDR));
    const kVals = rows.map(r => Number(r.k));
    const foldVals = rows.map(r => Number(r["Fold Enrichment"]));
    const geneSets = rows.map(_oraGeneSet);

    const geneIndex = {};
    const geneIds = [];
    geneSets.forEach((set) => {
        set.forEach((g) => {
            if (!(g in geneIndex)) { geneIndex[g] = "G" + geneIds.length; geneIds.push(g); }
        });
    });

    const nodeIds = termIds.concat(geneIds.map(g => geneIndex[g]));
    const edges = [];
    rows.forEach((r, i) => {
        geneSets[i].forEach((g) => edges.push({ source: termIds[i], target: geneIndex[g], weight: 1 }));
    });

    const pos = _forceLayout(nodeIds, edges, { iterations: 150 });

    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#D8E1EA", width: 1 }, hoverinfo: "skip", showlegend: false,
    };

    const geneLog2FCVals = geneIds.map(g => (g in geneLog2FC ? Number(geneLog2FC[g]) : null));
    const hasAnyFC = geneLog2FCVals.some(v => v !== null && isFinite(v));
    const maxAbsFC = hasAnyFC ? Math.max(0.5, ...geneLog2FCVals.filter(v => v !== null && isFinite(v)).map(Math.abs)) : 1;
    const geneTrace = {
        x: geneIds.map(g => pos[geneIndex[g]].x),
        y: geneIds.map(g => pos[geneIndex[g]].y),
        mode: "markers",
        type: "scatter",
        text: geneIds,
        customdata: geneLog2FCVals.map(v => (v === null || !isFinite(v)) ? "n/a" : v.toFixed(2)),
        hovertemplate: "%{text}<br>log2FC: %{customdata}<extra></extra>",
        marker: hasAnyFC ? {
            size: 9,
            color: geneLog2FCVals.map(v => (v === null || !isFinite(v)) ? 0 : v),
            colorscale: "RdBu", reversescale: true, cmid: 0,
            showscale: true,
            colorbar: { title: "log2FC", thickness: 12, len: 0.5, x: 1.12, y: 0.25 },
            line: { color: "#5D7086", width: 0.5 },
        } : { size: 8, color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } },
        showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const termTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], fdrVals[i], kVals[i], foldVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(kVals, 14, 38),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.5, y: 0.75 },
            line: { color: "#2c3e50", width: 1.5 },
            symbol: "diamond",
        },
        hovertemplate: "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Genes: %{customdata[2]}<br>Fold enrichment: %{customdata[3]:.2f}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, geneTrace, termTrace], {}, slug, null, "network", "Gene-concept network",
        () => _makeNetworkDraggable(
            div,
            [{ traceIndex: 1, nodeIds: geneIds.map(g => geneIndex[g]) }, { traceIndex: 2, nodeIds: termIds }],
            pos, 0, edges
        ));
}

// ---------------------------------------------------------------------------
// Burden Test relational charts: same two graphs as ORA (enrichment map,
// gene-concept network) but reading burden's columns ("Genes", "Effect (mean
// difference)", "Genes in set") and - in the gene-concept network - coloring
// each gene node by its variant score (from block.gene_scores). Reuse the same
// force layout / _plotNetwork / draggable plumbing as the ORA versions.
// ---------------------------------------------------------------------------

function buildBurdenEnrichmentMap(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const topN = opts.topN || 20;
    const fontSize = opts.fontSize || 12;
    const edgeThreshold = opts.edgeThreshold || 1;
    const rows = _burdenTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }
    if (rows.length < 2) { div.textContent = "Need at least 2 terms to draw an enrichment map."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const geneSets = rows.map(_burdenGeneSet);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r.FDR));
    const setSizes = rows.map(r => Number(r["Genes in set"]));
    const effectVals = rows.map(r => Number(r["Effect (mean difference)"]));

    const edges = [];
    for (let i = 0; i < rows.length; i++) {
        for (let j = i + 1; j < rows.length; j++) {
            let overlap = 0;
            geneSets[i].forEach((g) => { if (geneSets[j].has(g)) overlap++; });
            if (overlap < edgeThreshold) continue;
            const union = geneSets[i].size + geneSets[j].size - overlap;
            const jaccard = union > 0 ? overlap / union : 0;
            edges.push({ source: termIds[i], target: termIds[j], weight: jaccard, overlap });
        }
    }

    const pos = _forceLayout(termIds, edges, { iterations: 150 });

    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#B7C6D6", width: 1.5 }, hoverinfo: "skip", showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const nodeTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], fdrVals[i], setSizes[i], effectVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(setSizes, 10, 34),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.7 },
            line: { color: "#2c3e50", width: 1 },
        },
        hovertemplate: "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Genes in set: %{customdata[2]}<br>Effect: %{customdata[3]:.3f}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, nodeTrace], {}, dbName, namespace, "enrichmap", "Enrichment map",
        () => _makeNetworkDraggable(div, [{ traceIndex: 1, nodeIds: termIds }], pos, 0, edges));
}

function buildBurdenGeneConceptNetwork(div, table, dbName, namespace, colorScaleName, opts) {
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const geneScores = opts.geneScores || {};
    const topN = opts.topN === Infinity ? 15 : Math.min(opts.topN || 20, 15);
    const rows = _burdenTopRows(table, topN, opts.alpha);
    if (!rows.length) { div.textContent = "No significant terms to plot."; return; }

    const scale = ORA_DOT_COLOR_SCALES.find((item) => item[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
    const termIds = rows.map((r, i) => "T" + i);
    const names = rows.map(r => String(r["Term Name"] || ""));
    const fdrVals = rows.map(r => Number(r.FDR));
    const setSizes = rows.map(r => Number(r["Genes in set"]));
    const effectVals = rows.map(r => Number(r["Effect (mean difference)"]));
    const geneSets = rows.map(_burdenGeneSet);

    const geneIndex = {};
    const geneIds = [];
    geneSets.forEach((set) => {
        set.forEach((g) => {
            if (!(g in geneIndex)) { geneIndex[g] = "G" + geneIds.length; geneIds.push(g); }
        });
    });

    const nodeIds = termIds.concat(geneIds.map(g => geneIndex[g]));
    const edges = [];
    rows.forEach((r, i) => {
        geneSets[i].forEach((g) => edges.push({ source: termIds[i], target: geneIndex[g], weight: 1 }));
    });

    const pos = _forceLayout(nodeIds, edges, { iterations: 150 });

    const edgeX = [], edgeY = [];
    edges.forEach((e) => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = {
        x: edgeX, y: edgeY, mode: "lines", type: "scatter",
        line: { color: "#D8E1EA", width: 1 }, hoverinfo: "skip", showlegend: false,
    };

    // Gene nodes colored by their variant score (burden-specific). Genes with a
    // score get a second colorbar; any gene without one falls back to flat grey.
    const geneScoreVals = geneIds.map(g => (g in geneScores ? Number(geneScores[g]) : null));
    const hasAnyScore = geneScoreVals.some(v => v !== null && isFinite(v));
    const geneTrace = {
        x: geneIds.map(g => pos[geneIndex[g]].x),
        y: geneIds.map(g => pos[geneIndex[g]].y),
        mode: "markers",
        type: "scatter",
        text: geneIds,
        customdata: geneScoreVals.map(v => (v === null || !isFinite(v)) ? "n/a" : v),
        hovertemplate: "%{text}<br>Burden score: %{customdata}<extra></extra>",
        marker: hasAnyScore ? {
            size: 9,
            color: geneScoreVals.map(v => (v === null || !isFinite(v)) ? 0 : v),
            colorscale: "Viridis",
            showscale: true,
            colorbar: { title: "Burden score", thickness: 12, len: 0.5, x: 1.12, y: 0.25 },
            line: { color: "#5D7086", width: 0.5 },
        } : { size: 8, color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } },
        showlegend: false,
    };

    const wrappedNames = names.map(n => _wrapOraLabel(n, 28));
    const termTrace = {
        x: termIds.map(id => pos[id].x),
        y: termIds.map(id => pos[id].y),
        mode: "markers+text",
        type: "scatter",
        text: wrappedNames,
        textposition: "top center",
        textfont: { size: fontSize },
        customdata: rows.map((r, i) => [names[i], fdrVals[i], setSizes[i], effectVals[i]]),
        marker: {
            size: _magnitudeMarkerSizes(setSizes, 14, 38),
            color: fdrVals.map(f => -Math.log10(f + 1e-300)),
            colorscale: scale[2],
            reversescale: scale[3],
            showscale: true,
            colorbar: { title: "-log10(FDR)", thickness: 14, len: 0.5, y: 0.75 },
            line: { color: "#2c3e50", width: 1.5 },
            symbol: "diamond",
        },
        hovertemplate: "%{customdata[0]}<br>FDR: %{customdata[1]:.3e}<br>Genes in set: %{customdata[2]}<br>Effect: %{customdata[3]:.3f}<extra></extra>",
    };

    _plotNetwork(div, [edgeTrace, geneTrace, termTrace], {}, dbName, namespace, "network", "Gene-concept network",
        () => _makeNetworkDraggable(
            div,
            [{ traceIndex: 1, nodeIds: geneIds.map(g => geneIndex[g]) }, { traceIndex: 2, nodeIds: termIds }],
            pos, 0, edges
        ));
}

function renderOraResultBlock(dbName, namespace, block) {
    const dbCfg = STATE.ora.analysisConfig.perDatabase[dbName] || {};
    const alpha = Number(dbCfg.sigThreshold) > 0 ? Number(dbCfg.sigThreshold) : ORA_DEFAULT_SIG_THRESHOLD;

    const dbClass = dbName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
    const card = el("div", { class: "rank-card ora-db-card ora-result-card ora-db-card-" + dbClass });
    const title = namespace ? dbName + " — " + (ORA_GO_NAMESPACE_LABELS[namespace] || namespace) : dbName;
    card.appendChild(el("div", { class: "section-label ora-db-title" }, title));

    const navHost = el("div", {});
    card.appendChild(navHost);
    // stateKey includes dbName/namespace: each database gets its own card and
    // its own remembered section, rather than sharing one active-section slot.
    const nav = makeResultNav(navHost, { stateKey: "genomics_ora_" + dbName + "_" + (namespace || "x") });

    const overviewPane = nav.pane("Overview");
    if (block.caveat) overviewPane.appendChild(alertBox("warning", block.caveat));
    const cards = [
        ["Terms tested", block.n_terms_tested],
        ["Universe size", block.universe_size],
        ["Genes used", block.n_input_in_universe],
        ["Genes dropped", block.n_input_dropped],
    ];
    const recovered = block.n_recovered_by_symbol || 0;
    if (recovered > 0) cards.push(["Recovered by symbol", recovered]);
    overviewPane.appendChild(metricCards(cards, cards.length));
    if (recovered > 0) {
        overviewPane.appendChild(el("div", { class: "caption ora-symbol-recovery-note" },
            recovered + " gene(s) had no UniProt ID and were matched to GO by gene symbol (less precise than accession matching)."));
    }

    // Result table - its own pane, so no separate show/hide toggle is needed.
    const tablePane = nav.pane("Result table");
    tablePane.appendChild(buildTable(block.table, {
        rowClassFn: (row) => {
            const fdr = row["FDR"];
            return (fdr !== "" && fdr !== null && Number(fdr) < alpha) ? "row-significant" : null;
        },
        maxHeight: 480,
    }));
    const dlBtn = el("button", { class: "btn" }, "Download result table");
    dlBtn.addEventListener("click", () => {
        const path = "/api/ora/analysis/download/" + encodeURIComponent(dbName) + ".csv" + (namespace ? "?namespace=" + encodeURIComponent(namespace) : "");
        const filename = dbName + (namespace ? "_" + namespace : "") + "_ora_results.csv";
        downloadFile(path, filename);
    });
    tablePane.appendChild(dlBtn);

    const hasSigRows = (block.table.rows || []).some(row =>
        row["FDR"] !== null && row["FDR"] !== "" && isFinite(Number(row["FDR"])) && Number(row["FDR"]) < alpha
    );
    if (hasSigRows) {
        function makeChartPane(paneTitle) {
            const p = nav.pane(paneTitle, { group: "Charts" });
            const controls = el("div", { class: "chart-toggle-row" });
            const area = el("div", { class: "chart-pair" });
            p.appendChild(controls);
            p.appendChild(area);
            const toggles = [];
            return {
                el: p,
                option: (label, options, initialValue, onChange) => {
                    const control = el("label", { class: "chart-option-control" });
                    control.appendChild(el("span", { class: "chart-option-label" }, label));
                    const select = el("select", { class: "chart-option-select" });
                    options.forEach(([optionLabel, optionValue]) => {
                        select.appendChild(el("option", { value: optionValue }, optionLabel));
                    });
                    select.value = initialValue;
                    select.addEventListener("change", () => onChange(select.value));
                    control.appendChild(select);
                    controls.appendChild(control);
                    return select;
                },
                toggle: (kind, label, buildFn) => {
                    const chartId = _oraChartId(kind, dbName, namespace);
                    const isNetwork = kind === "enrichmap" || kind === "network";
                    const chartDiv = el("div", { id: chartId, class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
                    const chartBtn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
                    let chartShown = false;
                    let chartBuilt = false;
                    function renderChart() { buildFn(chartDiv); chartBuilt = true; }
                    chartBtn.addEventListener("click", () => {
                        if (!chartShown) {
                            chartDiv.style.display = "";
                            if (!chartBuilt) renderChart();
                            else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                            chartShown = true;
                            chartBtn.textContent = "Hide " + label;
                        } else {
                            chartDiv.style.display = "none";
                            chartShown = false;
                            chartBtn.textContent = "Show " + label;
                        }
                    });
                    controls.appendChild(chartBtn);
                    area.appendChild(chartDiv);
                    const t = { redraw: () => { if (chartShown) renderChart(); } };
                    toggles.push(t);
                    return t;
                },
                redraw: () => toggles.forEach(t => t.redraw()),
            };
        }

        // --- Gene set overview: dot, bar, volcano ---------------------------
        const overview = makeChartPane("Gene set overview");
        overview.el.insertBefore(el("div", { class: "ora-chart-threshold-note" },
            "Significance threshold: FDR < " + alpha + " (used for table highlighting and volcano-plot coloring; adjustable on the config screen)."),
            overview.el.firstChild);
        let dotColorScale1 = ORA_DOT_COLOR_SCALES[0][1];
        let barColor = ORA_BAR_COLORS[0][1];
        let topN = 20;
        let sortBy = "fdr";
        let fontSize = 12;
        let volcanoY = "p";  // volcano y-axis: "p" = -log10(p-value), "fdr" = -log10(FDR)

        // "All significant" maps to Infinity: _oraTopRows already gates every
        // chart to FDR < alpha first, so this shows every significant term
        // with no cap (rather than padding with non-significant ones).
        overview.option("Terms shown", [
            ["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"],
            ["50", "50"], ["100", "100"], ["200", "200"], ["All significant", "all"],
        ], String(topN === Infinity ? "all" : topN), (value) => {
            topN = value === "all" ? Infinity : Number(value);
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Rank & show by", [["FDR", "fdr"], ["Fold enrichment", "fold"], ["Gene count", "k"], ["Gene ratio", "ratio"]], sortBy, (value) => {
            sortBy = value;
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], String(fontSize), (value) => {
            fontSize = Number(value);
            overview.redraw();
            relationships.redraw();
        });

        const dotToggle = overview.toggle("dot", "dot plot", (div) => buildDotPlot(div, block.table, dbName, namespace, dotColorScale1, { topN, sortBy, fontSize, alpha }));
        overview.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale1, (value) => {
            dotColorScale1 = value;
            dotToggle.redraw();
        });
        const barToggle = overview.toggle("bar", "bar plot", (div) => buildBarPlot(div, block.table, dbName, namespace, barColor, { topN, sortBy, fontSize, alpha }));
        overview.option("Bar color", ORA_BAR_COLORS, barColor, (value) => {
            barColor = value;
            barToggle.redraw();
        });
        const volcanoToggle = overview.toggle("volcano", "volcano plot", (div) => buildVolcanoPlot(div, block.table, dbName, namespace, alpha, { fontSize, volcanoY }));
        overview.option("Volcano y-axis", [["-log10(p-value)", "p"], ["-log10(FDR)", "fdr"]], volcanoY, (value) => {
            volcanoY = value;
            volcanoToggle.redraw();
        });

        // --- Gene set relationships: enrichment map, gene-concept network ---
        const relationships = makeChartPane("Gene set relationships");
        relationships.el.insertBefore(el("div", { class: "ora-chart-threshold-note ora-chart-encoding-note" },
            "Circle color = -log10(FDR) (darker = more significant); circle size = number of your input genes in " +
            "the term. Rank & show by and Terms shown (Gene set overview) choose which terms are drawn."),
            relationships.el.firstChild);
        let dotColorScale2 = ORA_DOT_COLOR_SCALES[0][1];
        let edgeThreshold = 1;
        const enrichMapToggle = relationships.toggle("enrichmap", "enrichment map", (div) =>
            buildEnrichmentMap(div, block.table, dbName, namespace, dotColorScale2, { topN, sortBy, fontSize, edgeThreshold, alpha }));
        relationships.option("Edge threshold (shared genes)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (value) => {
            edgeThreshold = Number(value);
            enrichMapToggle.redraw();
        });
        const networkToggle = relationships.toggle("network", "gene-concept network", (div) =>
            buildGeneConceptNetwork(div, block.table, dbName, namespace, dotColorScale2, { topN, sortBy, fontSize, alpha }));
        relationships.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale2, (value) => {
            dotColorScale2 = value;
            enrichMapToggle.redraw();
            networkToggle.redraw();
        });
    }

    nav.finish();
    return card;
}

// ---------------------------------------------------------------------------
// Gene-set Burden Test results screen. Mirrors renderOraAnalysisResultsScreen
// / renderOraResultBlock structurally (metric cards, table, dot/bar/volcano
// charts with the same controls), reading burden's result columns and
// dispatching to the burden chart builders above. No enrichment-map/
// gene-concept-network here (out of scope for this analysis).
// ---------------------------------------------------------------------------

function renderBurdenAnalysisResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene-set Burden Test: results"));

    if (!STATE.ora.burdenResults) {
        wrap.appendChild(alertBox("info", "Running analysis..."));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "ora-analysis-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "ora-analysis-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    Object.keys(STATE.ora.burdenResults).forEach((dbName) => {
        const dbResult = STATE.ora.burdenResults[dbName];
        if (dbName === "GO") {
            ["biological_process", "molecular_function", "cellular_component"].forEach((ns) => {
                if (dbResult[ns]) wrap.appendChild(renderBurdenResultBlock(dbName, ns, dbResult[ns]));
            });
        } else {
            const key = Object.keys(dbResult)[0];
            if (key) wrap.appendChild(renderBurdenResultBlock(dbName, null, dbResult[key]));
        }
    });

    const newRunBtn = el("button", { class: "btn" }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.ora.burdenConfig = { databases: [], perDatabase: {} };
        STATE.ora.burdenResults = null;
        STATE.ora.burdenRunLocked = false;
        STATE.aaPage = "genomics_burden_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function renderBurdenResultBlock(dbName, namespace, block) {
    const dbCfg = STATE.ora.burdenConfig.perDatabase[dbName] || {};
    const alpha = Number(dbCfg.sigThreshold) > 0 ? Number(dbCfg.sigThreshold) : ORA_DEFAULT_SIG_THRESHOLD;

    const dbClass = dbName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
    const card = el("div", { class: "rank-card ora-db-card ora-result-card ora-db-card-" + dbClass });
    const title = namespace ? dbName + " — " + (ORA_GO_NAMESPACE_LABELS[namespace] || namespace) : dbName;
    card.appendChild(el("div", { class: "section-label ora-db-title" }, title));

    const navHost = el("div", {});
    card.appendChild(navHost);
    const nav = makeResultNav(navHost, { stateKey: "genomics_burden_" + dbName + "_" + (namespace || "x") });

    const overviewPane = nav.pane("Overview");
    if (block.caveat) overviewPane.appendChild(alertBox("warning", block.caveat));
    const cards = [
        ["Terms tested", block.n_terms_tested],
        ["Universe size", block.universe_size],
        ["Genes scored", block.n_input_in_universe],
        ["Genes dropped", block.n_input_dropped],
    ];
    overviewPane.appendChild(metricCards(cards, cards.length));

    // Orientation note: scores are always stored/tested so higher = more burden
    // (Effect > 0 and lower FDR both mean "more burden"). When the run used
    // "lower = more burden" (e.g. a raw p-value), the input was negated before
    // testing, so the Mean score columns below show negated values, not your
    // original numbers - flag this explicitly to avoid misreading them.
    const usedLowerDirection = STATE.ora.scoreSummary && STATE.ora.scoreSummary.direction === "lower";
    let orientationNote = "Effect > 0 and a lower FDR both mean the term's genes carry more burden.";
    if (usedLowerDirection) {
        orientationNote += " You selected \"lower = more burden\" for the score column, so the Mean score (set)/(rest) " +
            "columns below are shown negated from your original input values (still oriented so higher = more burden).";
    }
    overviewPane.appendChild(el("div", { class: "caption ora-burden-orientation-note" }, orientationNote));

    // Result table - its own pane, so no separate show/hide toggle is needed.
    const tablePane = nav.pane("Result table");
    tablePane.appendChild(buildTable(block.table, {
        rowClassFn: (row) => {
            const fdr = row["FDR"];
            return (fdr !== "" && fdr !== null && Number(fdr) < alpha) ? "row-significant" : null;
        },
        maxHeight: 480,
    }));
    const dlBtn = el("button", { class: "btn" }, "Download result table");
    dlBtn.addEventListener("click", () => {
        const path = "/api/burden/analysis/download/" + encodeURIComponent(dbName) + ".csv" + (namespace ? "?namespace=" + encodeURIComponent(namespace) : "");
        const filename = dbName + (namespace ? "_" + namespace : "") + "_burden_results.csv";
        downloadFile(path, filename);
    });
    tablePane.appendChild(dlBtn);

    const hasSigRows = (block.table.rows || []).some(row =>
        row["FDR"] !== null && row["FDR"] !== "" && isFinite(Number(row["FDR"])) && Number(row["FDR"]) < alpha
    );
    if (hasSigRows) {
        function makeChartPane(paneTitle) {
            const p = nav.pane(paneTitle, { group: "Charts" });
            const controls = el("div", { class: "chart-toggle-row" });
            const area = el("div", { class: "chart-pair" });
            p.appendChild(controls);
            p.appendChild(area);
            const toggles = [];
            return {
                el: p,
                option: (label, options, initialValue, onChange) => {
                    const control = el("label", { class: "chart-option-control" });
                    control.appendChild(el("span", { class: "chart-option-label" }, label));
                    const select = el("select", { class: "chart-option-select" });
                    options.forEach(([optionLabel, optionValue]) => {
                        select.appendChild(el("option", { value: optionValue }, optionLabel));
                    });
                    select.value = initialValue;
                    select.addEventListener("change", () => onChange(select.value));
                    control.appendChild(select);
                    controls.appendChild(control);
                    return select;
                },
                toggle: (kind, label, buildFn) => {
                    const chartId = _oraChartId(kind, dbName, namespace);
                    const isNetwork = kind === "enrichmap" || kind === "network";
                    const chartDiv = el("div", { id: chartId, class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
                    const chartBtn = el("button", { class: "btn chart-toggle-btn" }, "Show " + label);
                    let chartShown = false;
                    let chartBuilt = false;
                    function renderChart() { buildFn(chartDiv); chartBuilt = true; }
                    chartBtn.addEventListener("click", () => {
                        if (!chartShown) {
                            chartDiv.style.display = "";
                            if (!chartBuilt) renderChart();
                            else if (typeof Plotly !== "undefined") Plotly.Plots.resize(chartDiv);
                            chartShown = true;
                            chartBtn.textContent = "Hide " + label;
                        } else {
                            chartDiv.style.display = "none";
                            chartShown = false;
                            chartBtn.textContent = "Show " + label;
                        }
                    });
                    controls.appendChild(chartBtn);
                    area.appendChild(chartDiv);
                    const t = { redraw: () => { if (chartShown) renderChart(); } };
                    toggles.push(t);
                    return t;
                },
                redraw: () => toggles.forEach(t => t.redraw()),
            };
        }

        // --- Gene set overview: dot, bar, volcano ---------------------------
        const overview = makeChartPane("Gene set overview");
        overview.el.insertBefore(el("div", { class: "ora-chart-threshold-note" },
            "Significance threshold: FDR < " + alpha + " (used for table highlighting and volcano-plot coloring; adjustable on the config screen). " +
            "Dot x-axis = Effect (mean difference); bar length follows the Bar metric choice (FDR or Effect)."),
            overview.el.firstChild);
        let dotColorScale1 = ORA_DOT_COLOR_SCALES[0][1];
        let barColor = ORA_BAR_COLORS[0][1];
        let barMetric = "fdr";
        let topN = 20;
        let fontSize = 12;
        let volcanoY = "p";

        overview.option("Terms shown", [
            ["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"],
            ["50", "50"], ["100", "100"], ["200", "200"], ["All significant", "all"],
        ], String(topN === Infinity ? "all" : topN), (value) => {
            topN = value === "all" ? Infinity : Number(value);
            overview.redraw();
            relationships.redraw();
        });
        overview.option("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]], String(fontSize), (value) => {
            fontSize = Number(value);
            overview.redraw();
            relationships.redraw();
        });

        const dotToggle = overview.toggle("dot", "dot plot", (div) => buildBurdenDotPlot(div, block.table, dbName, namespace, dotColorScale1, { topN, fontSize, alpha }));
        overview.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale1, (value) => {
            dotColorScale1 = value;
            dotToggle.redraw();
        });
        const barToggle = overview.toggle("bar", "bar plot", (div) => buildBurdenBarPlot(div, block.table, dbName, namespace, barColor, { topN, fontSize, alpha, barMetric }));
        overview.option("Bar metric", [["FDR", "fdr"], ["Effect (mean difference)", "effect"]], barMetric, (value) => {
            barMetric = value;
            barToggle.redraw();
        });
        overview.option("Bar color", ORA_BAR_COLORS, barColor, (value) => {
            barColor = value;
            barToggle.redraw();
        });
        const volcanoToggle = overview.toggle("volcano", "volcano plot", (div) => buildBurdenVolcanoPlot(div, block.table, dbName, namespace, alpha, { fontSize, volcanoY }));
        overview.option("Volcano y-axis", [["-log10(p-value)", "p"], ["-log10(FDR)", "fdr"]], volcanoY, (value) => {
            volcanoY = value;
            volcanoToggle.redraw();
        });

        // --- Gene set relationships: enrichment map, gene-concept network ---
        // High-burden terms linked by shared genes, and which scored genes
        // drive them (gene nodes colored by score).
        const relationships = makeChartPane("Gene set relationships");
        relationships.el.insertBefore(el("div", { class: "ora-chart-threshold-note ora-chart-encoding-note" },
            "Term-circle color = -log10(FDR) and size = genes in set; in the gene-concept network the small gene " +
            "nodes are colored by their variant (burden) score."),
            relationships.el.firstChild);
        let dotColorScale2 = ORA_DOT_COLOR_SCALES[0][1];
        let edgeThreshold = 1;
        const enrichMapToggle = relationships.toggle("enrichmap", "enrichment map", (div) =>
            buildBurdenEnrichmentMap(div, block.table, dbName, namespace, dotColorScale2, { topN, fontSize, edgeThreshold, alpha }));
        relationships.option("Edge threshold (shared genes)", [[">= 1", "1"], [">= 2", "2"], [">= 3", "3"], [">= 5", "5"]], String(edgeThreshold), (value) => {
            edgeThreshold = Number(value);
            enrichMapToggle.redraw();
        });
        const networkToggle = relationships.toggle("network", "gene-concept network", (div) =>
            buildBurdenGeneConceptNetwork(div, block.table, dbName, namespace, dotColorScale2, { topN, fontSize, alpha, geneScores: block.gene_scores || {} }));
        relationships.option("Dot gradient", ORA_DOT_COLOR_SCALES.map(([label, value]) => [label, value]), dotColorScale2, (value) => {
            dotColorScale2 = value;
            enrichMapToggle.redraw();
            networkToggle.redraw();
        });
    }

    nav.finish();
    return card;
}

// ---------------------------------------------------------------------------
// Network Propagation / Module Detection: results screen (ADDITIVE by mode).
// The propagation block (summary cards + hit table + 2 figures) always shows;
// the module block (module table + module-colored network) is appended only
// when module detection was run. A "Parameters used" panel (with the fixed
// values too) is always shown for reproducibility.
// ---------------------------------------------------------------------------

// Distinct colors for module-membership coloring (recycled if >N modules).
const NETPROP_MODULE_COLORS = [
    "#1B639C", "#B84A52", "#2E8B57", "#7057A3", "#C77F1A", "#158A8A",
    "#B5539C", "#6C7A89", "#3D8B3D", "#A0522D", "#4682B4", "#8B4789",
];

// Selectable CATEGORICAL palettes for the module-colored network (one color per
// module). The propagation network uses a continuous gradient (ORA_DOT_COLOR_SCALES)
// instead, since it encodes a numeric score; a module is a category, so it needs
// a qualitative palette, not a gradient.
const NETPROP_MODULE_PALETTES = [
    ["Default", "default", NETPROP_MODULE_COLORS],
    ["Vivid", "vivid", ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]],
    ["Dark", "dark", ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]],
    ["Pastel", "pastel", ["#8dd3c7", "#b39ddb", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f"]],
];
function _netpropPalette(value) {
    const p = NETPROP_MODULE_PALETTES.find(x => x[1] === value);
    return p ? p[2] : NETPROP_MODULE_COLORS;
}

function _netpropSignificantHits(payload) {
    // The hit table carries every scored protein; the significant hits are
    // FDR < the run's threshold. Returns [{stringId, gene, score, z, fdr, seed}].
    const alpha = (payload.parameters && payload.parameters.fdr_alpha) || 0.05;
    const t = payload.propagation.table;
    return (t.rows || [])
        .filter(r => r["FDR"] !== null && isFinite(Number(r["FDR"])) && Number(r["FDR"]) < alpha)
        .map(r => ({
            stringId: r["STRING ID"], gene: r["Gene Symbol"] || r["STRING ID"],
            score: Number(r["Propagated score"]), z: Number(r["Z-score"]),
            fdr: Number(r["FDR"]), seed: r["Seed"] === true || r["Seed"] === "True",
        }));
}

function buildNetpropNetwork(div, payload, colorMode, opts) {
    // colorMode: "score" (heat by propagated score) or "module" (by cluster).
    // opts (all optional): nodesShown (top-N hits by propagated score; Infinity
    // = all), fontSize (node-label size), colorScaleName (score-mode gradient),
    // minEdgeConfidence (drop STRING links below this weight, 0-1), palette
    // (module-mode categorical color array).
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const nodesShown = opts.nodesShown || Infinity;
    const fontSize = opts.fontSize || 12;
    const colorScaleName = opts.colorScaleName || ORA_DOT_COLOR_SCALES[0][1];
    const minEdgeConfidence = opts.minEdgeConfidence || 0;
    const palette = (opts.palette && opts.palette.length) ? opts.palette : NETPROP_MODULE_COLORS;

    let hits = _netpropSignificantHits(payload);
    if (hits.length < 2) { div.textContent = "Not enough significant hits to draw a network."; return; }
    // Keep only the top-N hits by propagated score (the main declutter lever;
    // the full set is still in the hit table / CSV). Then the node set is those
    // top hits that share a drawn edge, so the layout stays clean.
    hits = hits.slice().sort((a, b) => b.score - a.score);
    if (isFinite(nodesShown)) hits = hits.slice(0, nodesShown);
    const edgesRaw = payload.propagation.edges || [];
    const info = {};
    hits.forEach(h => { info[h.stringId] = h; });

    const edgeKept = (e) => info[e.source] && info[e.target] && (e.weight == null || e.weight >= minEdgeConfidence);
    const nodeSet = new Set();
    edgesRaw.forEach(e => { if (edgeKept(e)) { nodeSet.add(e.source); nodeSet.add(e.target); } });
    if (nodeSet.size < 2) { div.textContent = "No STRING links to draw among these hits at the current settings (try more nodes or a lower minimum edge confidence)."; return; }
    const nodeIds = Array.from(nodeSet);
    const edges = edgesRaw
        .filter(edgeKept)
        .map(e => ({ source: e.source, target: e.target, weight: e.weight }));

    const pos = _forceLayout(nodeIds, edges, { iterations: 200 });

    const edgeX = [], edgeY = [];
    edges.forEach(e => { edgeX.push(pos[e.source].x, pos[e.target].x, null); edgeY.push(pos[e.source].y, pos[e.target].y, null); });
    const edgeTrace = { x: edgeX, y: edgeY, mode: "lines", type: "scatter", line: { color: "#C7D3E0", width: 1 }, hoverinfo: "skip", showlegend: false };

    const names = nodeIds.map(id => info[id].gene);
    const sizes = nodeIds.map(id => info[id].seed ? 15 : 10);
    const lineColors = nodeIds.map(id => info[id].seed ? "#111" : "#5D7086");
    const lineWidths = nodeIds.map(id => info[id].seed ? 2.2 : 0.7);

    let marker, hovertemplate, customdata;
    if (colorMode === "module") {
        const membership = (payload.modules && payload.modules.membership) || {};
        const modOf = {};
        Object.keys(membership).forEach((mid, i) => { (membership[mid] || []).forEach(sid => { modOf[sid] = { id: mid, idx: i }; }); });
        const colors = nodeIds.map(id => modOf[id] ? palette[modOf[id].idx % palette.length] : "#CBD5E0");
        customdata = nodeIds.map(id => [info[id].gene, modOf[id] ? modOf[id].id : "(none)", info[id].seed ? "seed" : "recruited"]);
        marker = { size: sizes, color: colors, line: { color: lineColors, width: lineWidths } };
        hovertemplate = "%{customdata[0]}<br>Module: %{customdata[1]}<br>%{customdata[2]}<extra></extra>";
    } else {
        const scaleName = ORA_DOT_COLOR_SCALES.find(s => s[1] === colorScaleName) || ORA_DOT_COLOR_SCALES[0];
        customdata = nodeIds.map(id => [info[id].gene, info[id].score, info[id].fdr, info[id].seed ? "seed" : "recruited"]);
        marker = {
            size: sizes,
            color: nodeIds.map(id => info[id].score),
            colorscale: scaleName[2], reversescale: scaleName[3], showscale: true,
            colorbar: { title: "Propagated<br>score", thickness: 14, len: 0.7 },
            line: { color: lineColors, width: lineWidths },
        };
        hovertemplate = "%{customdata[0]}<br>Propagated score: %{customdata[1]:.3e}<br>FDR: %{customdata[2]:.2e}<br>%{customdata[3]}<extra></extra>";
    }

    const nodeTrace = {
        x: nodeIds.map(id => pos[id].x), y: nodeIds.map(id => pos[id].y),
        mode: "markers+text", type: "scatter", text: names, textposition: "top center",
        textfont: { size: fontSize }, customdata, marker, hovertemplate, showlegend: false,
    };

    const kind = colorMode === "module" ? "modulenet" : "propnet";
    const label = colorMode === "module" ? "Module network" : "Propagation network";
    _plotNetwork(div, [edgeTrace, nodeTrace], {}, (opts.dbName || "netprop"), colorMode, kind, label,
        () => _makeNetworkDraggable(div, [{ traceIndex: 1, nodeIds }], pos, 0, edges));
}

function buildNetpropHistogram(div, payload, opts) {
    // Null-distribution view: histogram of every protein's Z-score, with the
    // ACTUAL significance boundary marked (the lowest Z-score among the
    // FDR-significant hits) - NOT a generic uncorrected z=1.645 line, which
    // would misrepresent the real (stricter, BH-FDR-based) cutoff used to
    // call hits. Shows WHY a protein counts as a hit (its signal is far above
    // the permutation null).
    if (typeof Plotly === "undefined") return;
    opts = opts || {};
    const fontSize = opts.fontSize || 12;
    const t = payload.propagation.table;
    const zs = (t.rows || []).map(r => Number(r["Z-score"])).filter(z => isFinite(z));
    if (!zs.length) { div.textContent = "No scores to plot."; return; }
    const alpha = (payload.parameters && payload.parameters.fdr_alpha) || 0.05;
    const hits = _netpropSignificantHits(payload);
    const trace = { x: zs, type: "histogram", nbinsx: 60, marker: { color: "#9DB2C6", line: { color: "#5D7086", width: 0.5 } }, hovertemplate: "Z-score %{x:.1f}<br>%{y} proteins<extra></extra>" };
    const layout = {
        font: { size: fontSize },
        xaxis: { title: { text: "Z-score (observed vs. permutation null)", font: { size: fontSize + 2 } }, tickfont: { size: fontSize } },
        yaxis: { title: { text: "Proteins", font: { size: fontSize + 2 } }, tickfont: { size: fontSize }, type: "log" },
        height: 420, margin: { l: 70, r: 30, t: 20, b: 56 },
        plot_bgcolor: "#FAFCFE", paper_bgcolor: "#FAFCFE",
    };
    if (hits.length) {
        const boundary = Math.min(...hits.map(h => h.z));
        layout.shapes = [{ type: "line", x0: boundary, x1: boundary, yref: "paper", y0: 0, y1: 1, line: { color: "#c0392b", width: 1.5, dash: "dash" } }];
        layout.annotations = [{ x: boundary, yref: "paper", y: 1, text: "FDR < " + alpha, showarrow: false, font: { size: Math.max(9, fontSize - 2), color: "#c0392b" }, xanchor: "left" }];
    }
    _newOraPlot(div, [trace], layout, (opts.dbName || "netprop"), "null", "histogram", "Null distribution");
}

function renderNetpropAnalysisResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Network Propagation / Module Detection: results"));

    const payload = STATE.ora.netpropResults;
    if (!payload) {
        wrap.appendChild(alertBox("info", "Running analysis... (loading the STRING network + permutation testing can take up to a minute or two)"));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "ora-analysis-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "ora-analysis-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    const prop = payload.propagation;
    const hasModules = !!payload.modules;

    // Bordered card with the standard blue top accent, matching GSEA/ORA -
    // see the identical comment in renderTranscriptomicsRegulonResultsScreen.
    const card = el("div", { class: "rank-card ora-db-card ora-result-card" });
    wrap.appendChild(card);
    const navHost = el("div", {});
    card.appendChild(navHost);
    const nav = makeResultNav(navHost, { stateKey: "genomics_netprop" });

    // ----- Always: propagation overview -----
    const overviewPane = nav.pane("Overview");
    overviewPane.appendChild(metricCards([
        ["Seeds mapped", prop.n_seeds_mapped],
        ["Seeds not in STRING", prop.n_seeds_not_in_string],
        ["Significant hits", prop.n_hits],
        ["Recruited (non-seed)", prop.n_recruited],
    ], 4));
    overviewPane.appendChild(el("div", { class: "caption" },
        "Network: " + prop.n_proteins + " proteins (" + (payload.taxon_name || "taxon " + payload.taxon_id) + ", STRING). " +
        "\"Recruited\" hits are proteins you did NOT upload that the network pulled in - proteins a standard enrichment of your uploaded gene list alone would never surface."));

    // ----- Hit table + download -----
    const hitPane = nav.pane("Hit table");
    hitPane.appendChild(buildTable(prop.table, {
        rowClassFn: (row) => (Number(row["FDR"]) < ((payload.parameters && payload.parameters.fdr_alpha) || 0.05) ? "row-significant" : null),
        maxHeight: 480,
        cellFormatFn: (col, val) => {
            if (col === "p-value" || col === "Empirical p-value" || col === "FDR"
                || col === "FDR (recruited family)") return _fmtSci(val, 3);
            if (col === "Propagated score") {
                if (val === "") return "";
                const num = Number(val);
                return isFinite(num) ? num.toExponential(2) : val;
            }
            if (col === "Z-score" || col === "Seed score") return _fmtFixed(val, 2);
            return val;
        },
    }));
    const hitDl = el("button", { class: "btn" }, "Download hit table (CSV)");
    hitDl.addEventListener("click", () => downloadFile("/api/netprop/analysis/download/hits.csv", "network_propagation_hits.csv"));
    hitPane.appendChild(hitDl);

    // Figures: propagation heat network + null histogram (always). All figures
    // share one controls row + state object (mirrors the ORA results
    // pattern): each control redraws whichever registered figures are shown.
    const figuresPane = nav.pane("Figures");
    const netpropFigState = { nodesShown: 50, fontSize: 12, colorScaleName: ORA_DOT_COLOR_SCALES[0][1], minEdgeConfidence: 0 };
    const netpropNetToggles = [];   // network figures: react to nodes / edge confidence
    const netpropAllToggles = [];   // all figures: react to font size
    let netpropScoreToggle = null;  // propagation (score) network: also reacts to color gradient
    const netpropControls = el("div", { class: "chart-toggle-row" });

    function netpropOption(label, options, initialValue, onChange) {
        const control = el("label", { class: "chart-option-control" });
        control.appendChild(el("span", { class: "chart-option-label" }, label));
        const select = el("select", { class: "chart-option-select" });
        options.forEach(([optLabel, optVal]) => select.appendChild(el("option", { value: optVal }, optLabel)));
        select.value = initialValue;
        select.addEventListener("change", () => onChange(select.value));
        control.appendChild(select);
        netpropControls.appendChild(control);
    }

    netpropOption("Nodes shown", [["20", "20"], ["30", "30"], ["50", "50"], ["100", "100"], ["200", "200"], ["All", "all"]],
        String(netpropFigState.nodesShown), (v) => { netpropFigState.nodesShown = v === "all" ? Infinity : Number(v); netpropNetToggles.forEach(t => t.redraw()); });
    netpropOption("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]],
        String(netpropFigState.fontSize), (v) => { netpropFigState.fontSize = Number(v); netpropAllToggles.forEach(t => t.redraw()); });
    netpropOption("Node color gradient", ORA_DOT_COLOR_SCALES.map(([lbl, val]) => [lbl, val]),
        netpropFigState.colorScaleName, (v) => { netpropFigState.colorScaleName = v; if (netpropScoreToggle) netpropScoreToggle.redraw(); });
    netpropOption("Min edge confidence", [["Any", "0"], [">= 0.7", "0.7"], [">= 0.8", "0.8"], [">= 0.9", "0.9"], [">= 0.95", "0.95"]],
        String(netpropFigState.minEdgeConfidence), (v) => { netpropFigState.minEdgeConfidence = Number(v); netpropNetToggles.forEach(t => t.redraw()); });
    figuresPane.appendChild(netpropControls);

    netpropScoreToggle = _netpropChartToggle("Show propagation network", "propnet", (div) => buildNetpropNetwork(div, payload, "score", netpropFigState));
    netpropNetToggles.push(netpropScoreToggle);
    netpropAllToggles.push(netpropScoreToggle);
    figuresPane.appendChild(netpropScoreToggle.el);

    const netpropHistToggle = _netpropChartToggle("Show null-distribution histogram", "histogram", (div) => buildNetpropHistogram(div, payload, netpropFigState));
    netpropAllToggles.push(netpropHistToggle);
    figuresPane.appendChild(netpropHistToggle.el);

    // ----- Conditional: module block -----
    if (hasModules) {
        const modulesPane = nav.pane("Modules");
        if (payload.modules.caveat) {
            // One banner for the whole section, not per module - every
            // module's enrichment shares the same universe, so a fetch/
            // propagation problem affects all of them identically.
            modulesPane.appendChild(alertBox("warning", payload.modules.caveat));
        }
        const mtable = payload.modules.table;
        modulesPane.appendChild(el("div", { class: "caption" }, (mtable.rows || []).length + " module(s) detected (MCL). Each is labeled by its top GO enrichment term."));
        modulesPane.appendChild(buildTable(mtable, {
            maxHeight: 400,
            cellFormatFn: (col, val) => {
                if (col === "Top enriched FDR" || col === "Top enriched FDR (raw, uncorrected)") return _fmtSci(val, 3);
                return val;
            },
        }));
        const modDl = el("button", { class: "btn" }, "Download module table (CSV)");
        modDl.addEventListener("click", () => downloadFile("/api/netprop/analysis/download/modules.csv", "network_modules.csv"));
        modulesPane.appendChild(modDl);
        // Module-colored network gets its OWN adjacent controls (independent of
        // the propagation network's row up top): font, nodes shown, min edge
        // confidence, and a categorical Module palette (its meaningful "color"
        // option, since it is colored by cluster rather than a numeric gradient).
        const moduleFigState = { nodesShown: 50, fontSize: 12, minEdgeConfidence: 0, palette: NETPROP_MODULE_PALETTES[0][1] };
        const moduleNetToggles = [];
        const moduleNetControls = el("div", { class: "chart-toggle-row" });
        function moduleNetOption(label, options, initialValue, onChange) {
            const control = el("label", { class: "chart-option-control" });
            control.appendChild(el("span", { class: "chart-option-label" }, label));
            const select = el("select", { class: "chart-option-select" });
            options.forEach(([optLabel, optVal]) => select.appendChild(el("option", { value: optVal }, optLabel)));
            select.value = initialValue;
            select.addEventListener("change", () => onChange(select.value));
            control.appendChild(select);
            moduleNetControls.appendChild(control);
        }
        moduleNetOption("Nodes shown", [["20", "20"], ["30", "30"], ["50", "50"], ["100", "100"], ["200", "200"], ["All", "all"]],
            String(moduleFigState.nodesShown), (v) => { moduleFigState.nodesShown = v === "all" ? Infinity : Number(v); moduleNetToggles.forEach(t => t.redraw()); });
        moduleNetOption("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]],
            String(moduleFigState.fontSize), (v) => { moduleFigState.fontSize = Number(v); moduleNetToggles.forEach(t => t.redraw()); });
        moduleNetOption("Module palette", NETPROP_MODULE_PALETTES.map(([lbl, val]) => [lbl, val]),
            moduleFigState.palette, (v) => { moduleFigState.palette = v; moduleNetToggles.forEach(t => t.redraw()); });
        moduleNetOption("Min edge confidence", [["Any", "0"], [">= 0.7", "0.7"], [">= 0.8", "0.8"], [">= 0.9", "0.9"], [">= 0.95", "0.95"]],
            String(moduleFigState.minEdgeConfidence), (v) => { moduleFigState.minEdgeConfidence = Number(v); moduleNetToggles.forEach(t => t.redraw()); });
        modulesPane.appendChild(moduleNetControls);

        const netpropModuleToggle = _netpropChartToggle("Show module-colored network", "modulenet",
            (div) => buildNetpropNetwork(div, payload, "module", {
                nodesShown: moduleFigState.nodesShown, fontSize: moduleFigState.fontSize,
                minEdgeConfidence: moduleFigState.minEdgeConfidence, palette: _netpropPalette(moduleFigState.palette),
            }));
        moduleNetToggles.push(netpropModuleToggle);
        modulesPane.appendChild(netpropModuleToggle.el);

        // Per-module functional enrichment: the module enrichment (ORA engine)
        // is computed in full, but the module table only shows each module's
        // single top term. Expose the rest as a dot/bar plot, one module at a
        // time (reuses the ORA dot/bar builders). Its own rail entry, since
        // it's a distinct interactive view (one module at a time), not part
        // of the Modules overview.
        const modEnrich = payload.modules.enrichment || {};
        const enrichModuleIds = (mtable.rows || [])
            .map(r => r["Module ID"])
            .filter(mid => modEnrich[mid] && modEnrich[mid].table && (modEnrich[mid].table.rows || []).length);
        if (enrichModuleIds.length) {
            const enrichPane = nav.pane("Per-module enrichment");
            enrichPane.appendChild(el("div", { class: "caption" },
                "Full functional enrichment for a chosen module (Fisher's exact + BH-FDR, the same engine as Over-Representation Analysis) - its top terms as a dot/bar plot, not just the single label in the table above. Circle/bar color = -log10(FDR); dot size = number of module genes in the term."));
            const modById = {};
            (mtable.rows || []).forEach(r => { modById[r["Module ID"]] = r; });
            let selectedEnrichId = enrichModuleIds[0];
            const enrichAlpha = (payload.parameters && payload.parameters.fdr_alpha) || 0.05;
            // Enrichment-plot state, matching the ORA enrichment controls. Font
            // size is its own control here (adjacent to the plots) rather than
            // the distant top-of-Figures one, so this section is self-contained.
            const enrichState = { topN: 15, sortBy: "fdr", fontSize: 12, dotColorScale: ORA_DOT_COLOR_SCALES[0][1], barColor: ORA_BAR_COLORS[0][1] };
            const enrichToggles = [];   // dot + bar; redraw on any enrichment control change

            const enrichControls = el("div", { class: "chart-toggle-row" });
            function enrichOption(label, options, initialValue, onChange) {
                const control = el("label", { class: "chart-option-control" });
                control.appendChild(el("span", { class: "chart-option-label" }, label));
                const select = el("select", { class: "chart-option-select" });
                options.forEach(([optLabel, optVal]) => select.appendChild(el("option", { value: optVal }, optLabel)));
                select.value = initialValue;
                select.addEventListener("change", () => onChange(select.value));
                control.appendChild(select);
                enrichControls.appendChild(control);
            }

            enrichOption("Module", enrichModuleIds.map(mid => {
                const r = modById[mid] || {};
                const term = r["Top enriched term"] ? " - " + r["Top enriched term"] : "";
                return [mid + " (" + (r["Size"] || "?") + " genes)" + term, mid];
            }), selectedEnrichId, (v) => { selectedEnrichId = v; enrichToggles.forEach(t => t.redraw()); });
            enrichOption("Terms shown", [["5", "5"], ["10", "10"], ["15", "15"], ["20", "20"], ["30", "30"], ["50", "50"], ["All", "all"]],
                String(enrichState.topN), (v) => { enrichState.topN = v === "all" ? Infinity : Number(v); enrichToggles.forEach(t => t.redraw()); });
            enrichOption("Rank & show by", [["FDR", "fdr"], ["Fold enrichment", "fold"], ["Gene count", "k"], ["Gene ratio", "ratio"]],
                enrichState.sortBy, (v) => { enrichState.sortBy = v; enrichToggles.forEach(t => t.redraw()); });
            enrichOption("Font size", [["10", "10"], ["12", "12"], ["14", "14"], ["16", "16"], ["18", "18"]],
                String(enrichState.fontSize), (v) => { enrichState.fontSize = Number(v); enrichToggles.forEach(t => t.redraw()); });
            enrichOption("Dot gradient", ORA_DOT_COLOR_SCALES.map(([lbl, val]) => [lbl, val]),
                enrichState.dotColorScale, (v) => { enrichState.dotColorScale = v; enrichToggles.forEach(t => t.redraw()); });
            enrichOption("Bar color", ORA_BAR_COLORS,
                enrichState.barColor, (v) => { enrichState.barColor = v; enrichToggles.forEach(t => t.redraw()); });
            enrichPane.appendChild(enrichControls);

            // dbName/namespace passed to buildDotPlot/buildBarPlot are STABLE
            // ("netprop_module"/"enrichment") - the container div is reused as
            // the selected module changes, so its id (and thus the export
            // lookup) must not depend on the module. The explicitId matches the
            // "dot"/"bar" kind those builders use internally for the export.
            const enrichDotToggle = _netpropChartToggle("Show enrichment dot plot", "moduleenrichdot",
                (div) => buildDotPlot(div, modEnrich[selectedEnrichId].table, "netprop_module", "enrichment",
                    enrichState.dotColorScale, { topN: enrichState.topN, sortBy: enrichState.sortBy, fontSize: enrichState.fontSize, alpha: enrichAlpha }),
                _oraChartId("dot", "netprop_module", "enrichment"));
            const enrichBarToggle = _netpropChartToggle("Show enrichment bar plot", "moduleenrichbar",
                (div) => buildBarPlot(div, modEnrich[selectedEnrichId].table, "netprop_module", "enrichment",
                    enrichState.barColor, { topN: enrichState.topN, sortBy: enrichState.sortBy, fontSize: enrichState.fontSize, alpha: enrichAlpha }),
                _oraChartId("bar", "netprop_module", "enrichment"));
            enrichToggles.push(enrichDotToggle, enrichBarToggle);
            enrichPane.appendChild(enrichDotToggle.el);
            enrichPane.appendChild(enrichBarToggle.el);
        }
    }

    // ----- Always: parameters used -----
    const paramsPane = nav.pane("Parameters used");
    const p = payload.parameters || {};
    const paramLines = [
        "STRING version: v" + (p.string_version || "12.0") + " (data licensed CC BY 4.0; cite STRING)",
        "STRING confidence cutoff: " + (p.confidence != null ? (p.confidence / 1000).toFixed(3) + " (combined_score >= " + p.confidence + ")" : "?"),
        "Restart probability (random walk): " + (p.restart_prob != null ? p.restart_prob : "0.5"),
        "Permutations: " + (p.n_permutations || "?") + " (" + (p.permutation_type || "degree-preserving") + ")",
        "Significance threshold: FDR < " + (p.fdr_alpha != null ? p.fdr_alpha : "0.05") + " (Benjamini-Hochberg on permutation z-score p-values)",
    ];
    if (hasModules) {
        paramLines.push("Module detection: Markov Clustering (MCL), inflation " + (p.mcl_inflation != null ? p.mcl_inflation : "2.0") + ", minimum module size " + (p.min_module_size || "3"));
        paramLines.push("Module enrichment: " + (p.enrichment_database || "GO") + " (Fisher's exact + BH-FDR, whole-genome background), term-size cap <=" + (p.module_max_term_size != null ? p.module_max_term_size : "2000") + " genes (excludes overly generic catch-all terms from module labeling)");
    }
    const paramBox = el("div", { class: "rank-card ora-db-card netprop-params-panel" });
    paramLines.forEach(line => paramBox.appendChild(el("div", { class: "caption" }, "- " + line)));
    paramsPane.appendChild(paramBox);

    nav.finish();

    const newRunBtn = el("button", { class: "btn" }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.ora.netpropResults = null;
        STATE.ora.netpropRunLocked = false;
        STATE.aaPage = "genomics_netprop_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

// Show/hide toggle for a netprop figure (builds lazily; reuses the ORA chart
// container styling + export via _plotNetwork/_newOraPlot inside buildFn).
// Returns { el, redraw }: redraw() rebuilds the figure (reading current shared
// control state via buildFn's closure) only if it is currently shown - mirrors
// ORA's addChartToggle so the controls row can re-render on option change.
function _netpropChartToggle(label, kind, buildFn, explicitId, dbName) {
    const holder = el("div", {});
    const isNetwork = kind === "propnet" || kind === "modulenet";
    const nsById = { histogram: "null", modulenet: "module", propnet: "score" };
    // explicitId lets a caller force the container id to match what its build
    // fn's export button will look up. The network/histogram builders pass a
    // `kind` that already matches; buildDotPlot/buildBarPlot hardcode kind
    // "dot"/"bar" internally, so the enrichment toggles MUST pass an explicitId
    // of _oraChartId("dot"|"bar", ...) or the modebar export can't find them.
    const chartDiv = el("div", { id: explicitId || _oraChartId(kind, dbName || "netprop", nsById[kind] || kind), class: "chart-container" + (isNetwork ? " chart-network" : ""), style: { display: "none" } });
    let shown = false, built = false;
    const btn = el("button", { class: "btn chart-toggle-btn" }, label);
    function renderChart() { buildFn(chartDiv); built = true; }
    btn.addEventListener("click", () => {
        if (!shown) {
            chartDiv.style.display = "";
            if (!built) { renderChart(); }
            else if (typeof Plotly !== "undefined") { Plotly.Plots.resize(chartDiv); }
            shown = true; btn.textContent = label.replace("Show", "Hide");
        } else { chartDiv.style.display = "none"; shown = false; btn.textContent = label; }
    });
    holder.appendChild(btn);
    holder.appendChild(chartDiv);
    return {
        el: holder,
        redraw: () => { if (shown) renderChart(); },
    };
}

// ---------------------------------------------------------------------------
// Gene Function Annotation: configure -> run -> results.
// No statistics, no significance threshold - the only setting is which GO
// evidence is shown by default (Include IEA), which the results screen can
// also flip LIVE (client-side, no re-run) since the engine always returns
// full evidence regardless of this setting (see Genomics_annotation.py).
// ---------------------------------------------------------------------------

function renderAnnotationAnalysisConfigScreen() {
    const wrap = el("div", { class: "ora-config-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene Function Annotation: setup"));

    const locked = STATE.ora.annotationRunLocked;
    if (locked) {
        wrap.appendChild(alertBox("info", "Setup is locked for the in-progress/completed run. Use \"Start new run\" on the results screen to change settings."));
    }

    wrap.appendChild(el("div", { class: "caption" },
        "This is a descriptive lookup, not a statistical test - there is no significance threshold to set. " +
        "The only choice is whether electronic (IEA) GO annotations are shown by default in the term-list " +
        "columns; you can flip this live on the results screen without re-running."));

    wrap.appendChild(el("div", { class: "section-label" }, "Parameters"));
    const paramBox = el("div", { class: "rank-card ora-db-card" });
    const row = el("div", { class: "ora-iea-row" });
    const cb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
    cb.checked = !!STATE.ora.annotationConfig.includeIea;
    cb.disabled = locked;
    cb.addEventListener("change", () => { STATE.ora.annotationConfig.includeIea = cb.checked; });
    row.appendChild(cb);
    row.appendChild(document.createTextNode(" Include IEA (electronic) GO annotations by default (adjustable live on the results screen)"));
    paramBox.appendChild(row);
    wrap.appendChild(paramBox);

    const actions = el("div", { class: "right-actions ora-config-actions" });
    const btnRow = el("div", { class: "btn-row ora-config-btn-row" });
    const backBtn = el("button", { class: "btn ora-back-btn" }, "Back");
    backBtn.addEventListener("click", () => { STATE.aaPage = "genomics_module"; render(); });
    btnRow.appendChild(backBtn);
    const runWrap = el("div", { class: "col" });
    const runBtn = el("button", { class: "btn btn-primary ora-run-btn" }, "Run");
    runBtn.disabled = locked || !STATE.ora.finalTable;
    runBtn.addEventListener("click", startAnnotationAnalysisRun);
    runWrap.appendChild(runBtn);
    btnRow.appendChild(runWrap);
    actions.appendChild(btnRow);
    wrap.appendChild(actions);

    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

function startAnnotationAnalysisRun() {
    if (!STATE.ora.finalTable) return;
    const payload = { include_iea: !!STATE.ora.annotationConfig.includeIea };
    apiPostJson("/api/annotation/analysis/configure", payload).then(() => {
        STATE.ora.annotationRunLocked = true;
        STATE.ora.annotationResults = null;
        STATE.aaPage = "genomics_annotation_results";
        render();
        return apiPostJson("/api/annotation/analysis/run", {});
    }).then(({ job_id }) => {
        pollJob(job_id, {
            onProgress: (job) => {
                const bar = document.getElementById("ora-analysis-progress-inner");
                const label = document.getElementById("ora-analysis-progress-label");
                if (bar) bar.style.width = Math.round((job.progress || 0) * 100) + "%";
                if (label) label.textContent = job.message || "";
            },
            onDone: (result) => { STATE.ora.annotationResults = result; render(); },
            onError: (err) => {
                alert("Gene function annotation run failed: " + err);
                STATE.ora.annotationRunLocked = false;
                render();
            },
        });
    }).catch((err) => {
        alert(err.message);
        STATE.ora.annotationRunLocked = false;
        render();
    });
}

const ANNOTATION_NS_TO_COLUMN = {
    "biological_process": "GO Biological Process",
    "molecular_function": "GO Molecular Function",
    "cellular_component": "GO Cellular Component",
};

function _annotationStatusClass(status) {
    if (status === "Well-characterized") return "annot-status-well";
    if (status === "Moderate") return "annot-status-moderate";
    if (status === "Sparse (dark)") return "annot-status-sparse";
    return null;
}

// Recomputes the collapsed table's GO term-list columns from the FULL-evidence
// tidy table when "Include IEA" is off - a pure client-side re-filter (no
// server round-trip), since the engine always returns full evidence. Every
// other column (incl. Annotation Status, Experimental/Total GO Terms) is left
// untouched, since those must always reflect full evidence regardless of the
// display toggle.
function _annotationFilteredCollapsedRows(payload, includeIea) {
    const baseRows = payload.collapsed.rows || [];
    if (includeIea) return baseRows;
    const byGene = {};
    (payload.tidy.rows || []).forEach(r => {
        if (r["Evidence Code"] === "IEA") return;
        const key = r["Gene"] + "||" + r["UniProt ID"];
        if (!byGene[key]) byGene[key] = { biological_process: new Set(), molecular_function: new Set(), cellular_component: new Set() };
        const bucket = byGene[key][r["Namespace"]];
        if (bucket && r["Term"]) bucket.add(r["Term"]);
    });
    return baseRows.map(row => {
        const key = row["Gene"] + "||" + row["UniProt ID"];
        const g = byGene[key] || { biological_process: new Set(), molecular_function: new Set(), cellular_component: new Set() };
        const newRow = Object.assign({}, row);
        Object.keys(ANNOTATION_NS_TO_COLUMN).forEach(ns => {
            newRow[ANNOTATION_NS_TO_COLUMN[ns]] = Array.from(g[ns]).sort().join("; ");
        });
        return newRow;
    });
}

// Splits one of the "; "-joined GO/Keywords cell strings back into its
// individual terms (empty/whitespace entries dropped).
function _annotationTermList(str) {
    return (str || "").split("; ").map(s => s.trim()).filter(Boolean);
}

function _annotationFieldRow(label, value) {
    return el("div", { class: "annot-field" }, [
        el("span", { class: "annot-field-label" }, label + ":"),
        el("span", { class: "annot-field-value" }, String(value)),
    ]);
}

// One click-to-open section inside a gene card (reuses the existing native
// details.expander component). opts.freeText renders `str` as a single
// paragraph with no count in the summary (used for UniProt Function); the
// default mode splits `str` on "; " and shows a term count.
function _annotationSection(label, str, opts) {
    opts = opts || {};
    const details = el("details", { class: "expander" });
    const body = el("div", { class: "expander-body" });
    if (opts.freeText) {
        const text = (str || "").trim();
        details.appendChild(el("summary", {}, label));
        body.appendChild(text ? el("div", {}, text) : el("div", { class: "annot-term-empty" }, "No curated function description."));
    } else {
        const items = _annotationTermList(str);
        details.appendChild(el("summary", {}, label + " (" + items.length + ")"));
        body.appendChild(items.length
            ? el("div", { class: "annot-term-list" }, items.join("; "))
            : el("div", { class: "annot-term-empty" }, "No terms."));
    }
    details.appendChild(body);
    return details;
}

// One per-gene "box": a labelled header (key-value fields, not table columns)
// plus five independent click-to-open sections for the bulky term lists /
// function text. `row` is one row from the (IEA-filtered) collapsed table.
function _annotationGeneCard(row) {
    const statusClass = _annotationStatusClass(row["Annotation Status"]) || "";
    const card = el("div", { class: "annot-gene-card " + statusClass });

    const header = el("div", { class: "annot-card-header" });
    header.appendChild(_annotationFieldRow("Gene", row["Gene"]));
    header.appendChild(_annotationFieldRow("UniProt ID", row["UniProt ID"]));
    header.appendChild(_annotationFieldRow("Protein Name", row["Protein Name"] || " - "));
    header.appendChild(_annotationFieldRow("Organism", row["Organism"] || " - "));
    header.appendChild(el("div", { class: "annot-field" }, [
        el("span", { class: "annot-field-label" }, "Annotation Status:"),
        el("span", { class: "annot-status-badge " + statusClass }, row["Annotation Status"]),
    ]));
    header.appendChild(_annotationFieldRow("Experimental GO Terms", row["Experimental GO Terms"]));
    header.appendChild(_annotationFieldRow("Total GO Terms", row["Total GO Terms"]));
    card.appendChild(header);

    const sections = el("div", { class: "annot-card-sections" });
    sections.appendChild(_annotationSection("GO Biological Process", row["GO Biological Process"]));
    sections.appendChild(_annotationSection("GO Molecular Function", row["GO Molecular Function"]));
    sections.appendChild(_annotationSection("GO Cellular Component", row["GO Cellular Component"]));
    sections.appendChild(_annotationSection("Keywords", row["Keywords"]));
    sections.appendChild(_annotationSection("UniProt Function", row["UniProt Function"], { freeText: true }));
    card.appendChild(sections);

    return card;
}

function renderAnnotationAnalysisResultsScreen() {
    const wrap = el("div", { class: "ora-config-screen ora-results-screen" });
    wrap.appendChild(aaHeader());
    wrap.appendChild(el("div", { class: "aa-omics-subtitle ora-config-title" }, "Genomics — Gene Function Annotation: results"));

    const payload = STATE.ora.annotationResults;
    if (!payload) {
        wrap.appendChild(alertBox("info", "Running analysis... (fetching GO + UniProt data can take a bit for larger gene lists)"));
        const progressWrap = el("div", { class: "progress-wrap" });
        const progressInner = el("div", { id: "ora-analysis-progress-inner", class: "progress-bar-inner", style: { width: "0%" } });
        progressWrap.appendChild(el("div", { class: "progress-bar-outer" }, progressInner));
        progressWrap.appendChild(el("div", { id: "ora-analysis-progress-label", class: "progress-label" }, ""));
        wrap.appendChild(progressWrap);
        return wrap;
    }

    // Bordered card with the standard blue top accent, matching GSEA/ORA -
    // see the identical comment in renderTranscriptomicsRegulonResultsScreen.
    const card = el("div", { class: "rank-card ora-db-card ora-result-card" });
    wrap.appendChild(card);
    const navHost = el("div", {});
    card.appendChild(navHost);
    const nav = makeResultNav(navHost, { stateKey: "genomics_annotation" });

    const counts = payload.counts || {};
    const overviewPane = nav.pane("Overview");
    overviewPane.appendChild(metricCards([
        ["Genes annotated", counts.n_genes || 0],
        ["Well-characterized", counts.well_characterized || 0],
        ["Moderate", counts.moderate || 0],
        ["Sparse (dark)", counts.sparse || 0],
    ], 4));
    overviewPane.appendChild(el("div", { class: "caption" },
        "\"Sparse (dark)\" genes have little experimental GO evidence and/or no curated UniProt function " +
        "description - they are the ones current databases know the least about, and are candidates for " +
        "further literature/experimental follow-up. See the README for the exact classification criteria."));

    const cardsPane = nav.pane("Gene cards");

    // Per-gene cards container, rebuilt LOCALLY on the IEA toggle (not a full
    // render()), so toggling doesn't jump scroll position or collapse other
    // open sections on the page.
    const cardsWrap = el("div", { class: "annot-gene-cards" });
    function rebuildCards() {
        const includeIea = !!STATE.ora.annotationConfig.includeIea;
        const displayRows = _annotationFilteredCollapsedRows(payload, includeIea);
        cardsWrap.innerHTML = "";
        displayRows.forEach((row) => cardsWrap.appendChild(_annotationGeneCard(row)));
    }

    // Live IEA toggle: re-filters the term-list sections instantly (no re-run).
    const toggleRow = el("div", { class: "ora-iea-row", style: { margin: "0.6rem 0" } });
    const toggleCb = el("input", { type: "checkbox", class: "ora-iea-checkbox" });
    toggleCb.checked = !!STATE.ora.annotationConfig.includeIea;
    toggleCb.addEventListener("change", () => {
        STATE.ora.annotationConfig.includeIea = toggleCb.checked;
        rebuildCards();
    });
    toggleRow.appendChild(toggleCb);
    toggleRow.appendChild(document.createTextNode(
        " Include IEA (electronic) GO annotations in the term lists below (instant - does not re-run the " +
        "lookup; Annotation Status always reflects full evidence either way)"
    ));
    cardsPane.appendChild(toggleRow);

    rebuildCards();
    cardsPane.appendChild(cardsWrap);

    const dlLinks = el("div", { class: "annot-download-links" });
    dlLinks.appendChild(document.createTextNode("Downloads: "));
    const dlAnnot = el("button", { class: "ora-inline-example-link", type: "button" }, "annotation table (CSV)");
    dlAnnot.addEventListener("click", () => downloadFile("/api/annotation/analysis/download/annotations.csv", "gene_function_annotations.csv"));
    const dlTidy = el("button", { class: "ora-inline-example-link", type: "button" }, "tidy term-level table (CSV)");
    dlTidy.addEventListener("click", () => downloadFile("/api/annotation/analysis/download/tidy.csv", "gene_function_annotations_tidy.csv"));
    dlLinks.appendChild(dlAnnot);
    dlLinks.appendChild(document.createTextNode(" | "));
    dlLinks.appendChild(dlTidy);
    cardsPane.appendChild(dlLinks);

    nav.finish();

    const newRunBtn = el("button", { class: "btn", style: { marginTop: "0.8rem" } }, "Start new run");
    newRunBtn.addEventListener("click", () => {
        STATE.ora.annotationResults = null;
        STATE.ora.annotationRunLocked = false;
        STATE.aaPage = "genomics_annotation_config";
        render();
    });
    wrap.appendChild(newRunBtn);
    wrap.appendChild(el("div", { class: "aa-bottom-nav-spacer" }));
    wrap.appendChild(backToGenomicsActions());
    return wrap;
}

document.addEventListener("DOMContentLoaded", () => {
    STATE.ora = freshOraState();
    STATE.transcriptomics = freshTranscriptomicsState();
    render();
});
