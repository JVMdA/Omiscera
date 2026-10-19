# Omiscera

Omiscera is a web application for **annotation analytics across omics layers**: Genomics, Transcriptomics,
Proteomics, Metabolomics, and six Multi-omics integration methods, all in one tool. It is built for researchers who
are not necessarily bioinformaticians — every result explains the statistics behind it, states its assumptions, and
says plainly what it does *not* mean.

A FastAPI backend (`API.py` + one module per analysis, e.g. `Transcriptomics_module.py`, `Multiomics_mofa.py`) serves
a vanilla-JavaScript single-page frontend (`static/app.js`). There is no build step for the frontend.

## What it does

- **Genomics** — over-representation analysis (ORA), network propagation, burden testing, gene function annotation.
- **Transcriptomics** — ORA, GSEA, network/module detection, KEGG mapping, TF–regulon analysis.
- **Proteomics** — ORA, quantitative enrichment, PPI network, PTM analysis.
- **Metabolomics** — ORA, quantitative enrichment, pathway topology analysis.
- **Multi-omics** — Joint Pathway Enrichment, Cross-omics Pathway Mapping, Correlation-based Integration,
  Network-based Integration, MOFA (Factor-based Integration), and DIABLO (Multiblock Integration).

Each module has its own README under [`info/`](info/) explaining the method in plain terms, its statistics, its
inputs, and what its results do and don't tell you — read the one for the method you're using before interpreting
its output.

## Requirements

- **Python 3.12.** No other language or runtime is needed to run Omiscera (R is only ever used to regenerate one
  test fixture for the DIABLO method — see [Tests](#tests) — never to run the app itself).
- Every Python package Omiscera needs, pinned to an exact version, is in [`requirements.txt`](requirements.txt):
  `numpy`, `scipy`, `pandas`, `scikit-learn`, `h5py` and `mofapy2` for the statistics and the MOFA method;
  `fastapi`, `starlette`, `pydantic`, `uvicorn`, `python-multipart` and `anyio` for the web server; `requests`,
  `urllib3`, `certifi`, `idna` for talking to the external databases listed on the app's own References page;
  `matplotlib`, `pillow` and `openpyxl` for figure export and spreadsheet input/output. `playwright` is listed too
  but is only needed to run the browser-driven tests, not to run the app.
- No database server, no separate frontend build step, and no extra system libraries to install — `pip install -r
  requirements.txt` is everything.
- Works on Windows, macOS and Linux. The one platform-specific gotcha: `numpy` must come from a standard pip build
  (see below) — this was found to matter on Windows with a bundled scientific Python distribution, not on a plain
  `python.org`/`pyenv` install.

## Quickstart

Create a virtual environment and install the pinned dependencies:

```bash
python3.12 -m venv .venv
# Windows:
.venv\Scripts\python.exe -m pip install -r requirements.txt
# macOS / Linux:
.venv/bin/python -m pip install -r requirements.txt
```

`numpy` must come from a standard pip build (it ships its own BLAS) — some bundled scientific Python distributions
have numpy installs with no working BLAS, which makes every matrix operation crash; MOFA and DIABLO both need a
working BLAS. If a bundled distribution's Python gives you crashes on basic linear algebra, use a plain `python.org`
or `pyenv` install instead.

Start the server:

```bash
.venv/bin/python -m uvicorn API:app --port 8000
# or just:
.venv/bin/python API.py
```

Open `http://127.0.0.1:8000` in a browser.

Browser-driven tests additionally need Playwright's Chromium:

```bash
.venv/bin/python -m playwright install chromium
```

## Data and caching

External reference data (listed on the app's own References page — STRING, Reactome, GO, MSigDB, GWAS Catalog, ...)
is fetched on demand and cached under `cache/` (not tracked in git — about 1.2 GB once fully populated, and fully
regenerable). `Data2/` and `Data3/` hold small, bundled example input files the app's "Example file" links and the
test suite use.

## Tests

The suite is plain Python scripts, not pytest — each one does its own pass/fail bookkeeping and exits 0 or 1:

```bash
python tests/test_regulon_engine.py
```

See [`tests/README.txt`](tests/README.txt) for what each file covers, which need a running server (port 8009 by
convention for that suite) or network access, and how the mixOmics-agreement fixtures for the DIABLO engine are
regenerated (`tests/r_reference/README.txt` — only needed if you touch `Multiomics_diablo.py`, not to run Omiscera).

## Citation

If you use Omiscera in your research, please cite it — see [`CITATION.cff`](CITATION.cff) (GitHub turns this into a
"Cite this repository" button automatically).

## License

[MIT](LICENSE).
