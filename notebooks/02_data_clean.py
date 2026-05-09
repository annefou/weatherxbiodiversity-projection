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
# # 02 — Data clean (Pass 1, Iberian baseline)
#
# This notebook is a **thin orchestration wrapper** around the vendored
# upstream pipeline in `soroye_port/`. It runs four upstream scripts in
# order:
#
# 1. `01_clean_data_iberia.py` — clean GBIF Iberia download, apply
#    Kerr-2015 species filter and IUCN exclusion list.
# 2. `02_presence_absence.py` — build the 100km cylindrical-equal-area
#    grid spanning N. America + Europe and infer presence/absence per
#    (species × period × season).
# 3. `03_sampling_continent.py` — per-cell sampling-effort raster
#    (distinct LYIDs) and the continent code.
# 4. `04_climate_tei_pei.py` — bilinearly interpolate CRU TS 3.24.01
#    onto the CEA grid and compute Thermal & Precipitation Exposure
#    Indices per species per cell per period.
#
# We do **not** re-port the science. The upstream scripts produced
# `sc_TEI_delta = +0.479 [0.265, 0.694]` in v0.2.1 and we want to
# reproduce that headline byte-for-byte before changing anything. The
# only environment toggle is `OUT_SUBDIR=outputs_iberia`, which selects
# the Phase-3 (Iberia) output directory inside `soroye_port/`.

# %%
import os
import subprocess
import sys
from pathlib import Path

# %%
ROOT = Path("..").resolve()
PORT = ROOT / "soroye_port"
OUT_DIR = PORT / "outputs_iberia"

env = {**os.environ, "OUT_SUBDIR": "outputs_iberia"}


def run(script: str) -> None:
    print(f"\n=== {script} ===", flush=True)
    subprocess.run(
        [sys.executable, script],
        cwd=PORT,
        env=env,
        check=True,
    )


# %% [markdown]
# ## Run the four upstream cleaning + indexing scripts

# %%
run("01_clean_data_iberia.py")
run("02_presence_absence.py")
run("03_sampling_continent.py")
run("04_climate_tei_pei.py")

# %% [markdown]
# ## Summary of intermediate artefacts produced

# %%
expected = [
    OUT_DIR / "bombus_clean.csv",
    OUT_DIR / "presence_absence.npz",
    OUT_DIR / "sampling_continent.npz",
    OUT_DIR / "climate_tei_pei.npz",
]
print("\nIntermediate artefacts:")
for p in expected:
    if p.exists():
        size = p.stat().st_size
        print(f"  ok    {p.relative_to(ROOT)}  ({size:,} bytes)")
    else:
        print(f"  MISS  {p.relative_to(ROOT)}")

# %%
# Quick row-count check on the cleaned occurrence table — gives the
# user a one-line confirmation that the pipeline reached a sensible
# place before the regression runs in 03_analysis.py.
import pandas as pd  # noqa: E402

clean_csv = OUT_DIR / "bombus_clean.csv"
if clean_csv.exists():
    df = pd.read_csv(clean_csv)
    print(f"\nbombus_clean.csv → {len(df):,} rows, {df['species'].nunique()} species")
    print(f"  periods present: {sorted(df['period'].dropna().unique().tolist())}")
    print(f"  seasons present: {sorted(df['season'].dropna().unique().tolist())}")
