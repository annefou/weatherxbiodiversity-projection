# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 01 — Data download (Pass 1, Iberian baseline)
#
# This notebook fetches **all** input data needed by the vendored
# upstream pipeline in `soroye_port/`. There are three artefacts:
#
# 1. **GBIF Iberia *Bombus* download** — fetched directly from the
#    pre-existing GBIF download key `0006204-260423192947929` (DOI
#    [`10.15468/dl.3frmsq`](https://doi.org/10.15468/dl.3frmsq)). We do
#    not re-mint a new download — this is the exact snapshot the
#    upstream's `posterior_vb_summary.csv` was built from, and re-minting
#    would introduce a fresh DOI and change the GBIF citation
#    provenance. No GBIF credentials needed.
# 2. **Soroye et al. Figshare deposit** ([10.6084/m9.figshare.9956471](https://doi.org/10.6084/m9.figshare.9956471))
#    — a ~1.13 GB ZIP containing the Kerr 2015 species list, the
#    erroneous-observations table, and the bundled CRU TS 3.24.01
#    monthly climate NetCDFs. The vendored
#    `soroye_port/04_climate_tei_pei.py` reads these NetCDFs directly.
# 3. **GBIF citation metadata** — a small JSON written into `data/`
#    recording the download key, DOI and citation string for downstream
#    provenance.
#
# All three are cached: re-running the notebook only re-downloads if a
# file is missing or truncated.

# %%
import json
import os
import shutil
import zipfile
from pathlib import Path

import requests

# %%
# Paths — match exactly what the vendored upstream scripts hard-code.
ROOT = Path("..").resolve()
DATA_DIR = ROOT / "data"
GBIF_DL_DIR = DATA_DIR / "gbif_dl"
REF_DIR = ROOT / "reference"
GBIF_DL_DIR.mkdir(parents=True, exist_ok=True)
REF_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## 1. Pre-minted GBIF download (DOI `10.15468/dl.3frmsq`)
#
# The upstream `01_clean_data_iberia.py` reads a TSV at
# `data/gbif_dl/0006204-260423192947929.csv` (despite the `.csv`
# extension, GBIF's SIMPLE_CSV format is tab-separated). We fetch the
# zipped download from the public GBIF API endpoint and extract the
# inner TSV under that exact filename.

# %%
GBIF_DL_KEY = "0006204-260423192947929"
GBIF_DL_URL = f"https://api.gbif.org/v1/occurrence/download/request/{GBIF_DL_KEY}.zip"
GBIF_DL_DOI = "10.15468/dl.3frmsq"

GBIF_ZIP_PATH = GBIF_DL_DIR / f"{GBIF_DL_KEY}.zip"
GBIF_TSV_PATH = GBIF_DL_DIR / f"{GBIF_DL_KEY}.csv"

# Conservative size threshold: the zip is ~1.9 MB, so >100 KB is "looks real".
MIN_ZIP_BYTES = 100_000
MIN_TSV_BYTES = 1_000_000


def fetch_gbif_zip() -> Path:
    if GBIF_ZIP_PATH.exists() and GBIF_ZIP_PATH.stat().st_size > MIN_ZIP_BYTES:
        print(f"  cached: {GBIF_ZIP_PATH} ({GBIF_ZIP_PATH.stat().st_size:,} bytes)")
        return GBIF_ZIP_PATH
    print(f"  fetching {GBIF_DL_URL}")
    r = requests.get(GBIF_DL_URL, stream=True, timeout=600, allow_redirects=True)
    r.raise_for_status()
    with open(GBIF_ZIP_PATH, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 20):
            f.write(chunk)
    print(f"  saved {GBIF_ZIP_PATH} ({GBIF_ZIP_PATH.stat().st_size:,} bytes)")
    return GBIF_ZIP_PATH


def extract_gbif_tsv() -> Path:
    if GBIF_TSV_PATH.exists() and GBIF_TSV_PATH.stat().st_size > MIN_TSV_BYTES:
        print(f"  cached: {GBIF_TSV_PATH} ({GBIF_TSV_PATH.stat().st_size:,} bytes)")
        return GBIF_TSV_PATH
    with zipfile.ZipFile(GBIF_ZIP_PATH) as zf:
        members = zf.namelist()
        # SIMPLE_CSV downloads contain one occurrence file (.csv extension,
        # but tab-separated). Find it regardless of internal name.
        candidates = [m for m in members if m.endswith(".csv") or m.endswith(".tsv")]
        if not candidates:
            raise RuntimeError(
                f"No CSV/TSV inside {GBIF_ZIP_PATH}; members={members}"
            )
        # Pick the largest (avoids picking a manifest)
        candidates.sort(key=lambda m: zf.getinfo(m).file_size, reverse=True)
        member = candidates[0]
        print(f"  extracting {member} → {GBIF_TSV_PATH}")
        with zf.open(member) as src, open(GBIF_TSV_PATH, "wb") as dst:
            shutil.copyfileobj(src, dst)
    print(f"  saved {GBIF_TSV_PATH} ({GBIF_TSV_PATH.stat().st_size:,} bytes)")
    return GBIF_TSV_PATH


print("--- 1. GBIF Iberia Bombus download ---")
fetch_gbif_zip()
extract_gbif_tsv()

# %% [markdown]
# ## 2. Soroye Figshare deposit (Kerr species list + CRU TS climate)
#
# The figshare API resolves
# [10.6084/m9.figshare.9956471](https://doi.org/10.6084/m9.figshare.9956471)
# to a single `Bumblebee_repo.zip` (~1.13 GB). Extracted, it produces
# `Bumblebee_repo_wbombusdat/` containing both the species filter
# tables and the CRU TS 3.24.01 monthly NetCDFs that
# `soroye_port/04_climate_tei_pei.py` opens directly.

# %%
FIGSHARE_API = "https://api.figshare.com/v2/articles/9956471"
BUMBLEBEE_ZIP = REF_DIR / "Bumblebee_repo.zip"
BUMBLEBEE_DIR = REF_DIR / "Bumblebee_repo_wbombusdat"

# Files the vendored upstream scripts expect to exist after extraction.
EXPECTED_FILES = [
    BUMBLEBEE_DIR / "0_data" / "Kerr_et_al2015_spplist.csv",
    BUMBLEBEE_DIR / "0_data" / "bombus_err_obs.csv",
]

# CRU TS climate dir — `04_climate_tei_pei.py` globs `*tmp.dat.nc`,
# `*tmn.dat.nc`, `*tmx.dat.nc`, `*pre.dat.nc` from this directory.
CRU_DIR = BUMBLEBEE_DIR / "0_ClimateData"

# Conservative threshold for the zip: > 100 MB means "actually large".
MIN_BUMBLEBEE_ZIP_BYTES = 100_000_000


def figshare_download_url() -> str:
    print(f"  querying {FIGSHARE_API}")
    r = requests.get(FIGSHARE_API, timeout=120)
    r.raise_for_status()
    meta = r.json()
    files = meta.get("files", [])
    if not files:
        raise RuntimeError(f"No files in figshare article: {meta}")
    # The deposit currently has a single Bumblebee_repo.zip.
    target = next(
        (f for f in files if f["name"].lower().startswith("bumblebee")),
        files[0],
    )
    print(f"  resolved {target['name']} ({target['size']:,} bytes)")
    return target["download_url"]


def fetch_bumblebee_zip() -> Path:
    if BUMBLEBEE_ZIP.exists() and BUMBLEBEE_ZIP.stat().st_size > MIN_BUMBLEBEE_ZIP_BYTES:
        print(f"  cached: {BUMBLEBEE_ZIP} ({BUMBLEBEE_ZIP.stat().st_size:,} bytes)")
        return BUMBLEBEE_ZIP
    url = figshare_download_url()
    print(f"  fetching {url} (this is large, ~1.1 GB)")
    r = requests.get(url, stream=True, timeout=3600, allow_redirects=True)
    r.raise_for_status()
    written = 0
    with open(BUMBLEBEE_ZIP, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 20):
            f.write(chunk)
            written += len(chunk)
            if written % (50 << 20) < (1 << 20):  # progress every ~50 MB
                print(f"    {written / (1 << 20):,.0f} MB")
    print(f"  saved {BUMBLEBEE_ZIP} ({BUMBLEBEE_ZIP.stat().st_size:,} bytes)")
    return BUMBLEBEE_ZIP


def extract_bumblebee_zip() -> Path:
    if all(p.exists() for p in EXPECTED_FILES) and CRU_DIR.exists():
        print(f"  cached: extracted tree under {BUMBLEBEE_DIR}")
        return BUMBLEBEE_DIR
    print(f"  extracting {BUMBLEBEE_ZIP} → {REF_DIR}")
    with zipfile.ZipFile(BUMBLEBEE_ZIP) as zf:
        zf.extractall(REF_DIR)
    # Validate
    missing = [p for p in EXPECTED_FILES if not p.exists()]
    if missing:
        raise RuntimeError(
            f"Extraction completed but expected files are missing: {missing}. "
            f"Inspect {BUMBLEBEE_DIR} to discover the correct sub-paths."
        )
    if not CRU_DIR.exists():
        raise RuntimeError(
            f"Expected CRU TS climate directory not found at {CRU_DIR}. "
            f"04_climate_tei_pei.py needs *.tmp.dat.nc / *.tmn.dat.nc / "
            f"*.tmx.dat.nc / *.pre.dat.nc files in this directory."
        )
    print(f"  validated expected files under {BUMBLEBEE_DIR}")
    return BUMBLEBEE_DIR


print("\n--- 2. Soroye Figshare deposit (Kerr species list + CRU TS) ---")
fetch_bumblebee_zip()
extract_bumblebee_zip()

# Also report the discovered CRU TS NetCDFs so the user sees what
# 04_climate_tei_pei.py will load.
if CRU_DIR.exists():
    cru_files = sorted(CRU_DIR.glob("*.dat.nc"))
    print(f"  CRU TS NetCDFs in {CRU_DIR.name}: {len(cru_files)}")
    for p in cru_files[:8]:
        print(f"    {p.name}")
    if len(cru_files) > 8:
        print(f"    … and {len(cru_files) - 8} more")

# %% [markdown]
# ## 3. GBIF citation metadata
#
# A small JSON record so downstream provenance (and the FORRT
# Replication Study draft) can quote the GBIF citation verbatim.
# This file is intentionally tracked in git (`data/.gitignore` has an
# exception) — it is the citation contract for this repository.

# %%
GBIF_META_PATH = DATA_DIR / "gbif_bombus_iberia_metadata.json"

gbif_meta = {
    "download_key": GBIF_DL_KEY,
    "doi": GBIF_DL_DOI,
    "doi_url": f"https://doi.org/{GBIF_DL_DOI}",
    "citation": (
        "GBIF.org (2026-04-25) GBIF Occurrence Download "
        f"https://doi.org/{GBIF_DL_DOI}"
    ),
    "source_url": GBIF_DL_URL,
    "format": "SIMPLE_CSV (tab-separated)",
    "local_path": str(GBIF_TSV_PATH.relative_to(ROOT)),
    "license": "CC-BY-NC-4.0 (per individual GBIF datasets)",
    "upstream_repo": {
        "name": "annefou/weatherxbiodiversity",
        "version": "v0.2.1",
        "concept_doi": "10.5281/zenodo.19756173",
    },
}
with open(GBIF_META_PATH, "w") as f:
    json.dump(gbif_meta, f, indent=2)
print(f"\n--- 3. Wrote {GBIF_META_PATH}")
print(f"    citation: {gbif_meta['citation']}")

# %% [markdown]
# ## Summary
#
# All three downloads are now in place. The vendored upstream pipeline
# in `soroye_port/` will read them at the exact paths it hard-codes.

# %%
artefacts = [
    ("GBIF zip",     GBIF_ZIP_PATH),
    ("GBIF TSV",     GBIF_TSV_PATH),
    ("Bumblebee zip", BUMBLEBEE_ZIP),
    ("Kerr species list", EXPECTED_FILES[0]),
    ("err-obs table", EXPECTED_FILES[1]),
    ("CRU TS dir",   CRU_DIR),
    ("citation JSON", GBIF_META_PATH),
]
for name, p in artefacts:
    if p.exists():
        size = p.stat().st_size if p.is_file() else sum(
            f.stat().st_size for f in p.rglob("*") if f.is_file()
        )
        print(f"  ok    {name:<22} {size:>15,} bytes  {p.relative_to(ROOT)}")
    else:
        print(f"  MISS  {name:<22} {'?':>15}        {p.relative_to(ROOT)}")
