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

# %%
import _tier2_guard; _tier2_guard.ensure_destine_or_skip()

# %% [markdown]
# # 05 — DestinE Climate DT extract (Tier 2)
#
# Fetches **SSP3-7.0** daily fields (max temperature, min temperature,
# total precipitation) over the Iberian Peninsula bounding box for two
# decade slices:
#
#   * **Mid-century**: 2046–2055
#   * **End-of-century**: 2076–2085
#
# Output NetCDFs land in `data/destine/` and are **gitignored** — the
# DestinE Climate DT data is licence-locked (no redistribution); only
# aggregated derived statistics may be committed.
#
# This notebook only runs on the DestinE Jupyter platform (or any
# environment with valid `polytope-client` credentials at
# `~/.polytopeapirc`). On every other host the import-time guard exits
# the notebook cleanly.
#
# **CHECK ITEMS** for first-run on the DestinE platform:
#   1. Verify the catalogue keys (class / dataset / activity / experiment
#      / model / resolution / type / stream / param) below match the
#      current DestinE Climate DT collection. Param numbers `165/166/167`
#      are placeholders for max temperature / min temperature / total
#      precipitation in the historical-era ECMWF GRIB-1 table — the
#      actual DestinE Climate DT param identifiers may differ.
#   2. Confirm the time-step granularity (daily vs hourly) and adjust
#      the date range syntax accordingly.
#   3. Confirm the Iberia bounding-box clipping happens server-side
#      (``area``) or client-side (``ds.sel(lat=..., lon=...)``).

# %%
import os
from pathlib import Path

import xarray as xr

# %%
ROOT = Path("..").resolve()
DATA_DIR = ROOT / "data" / "destine"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Iberian Peninsula bounding box (deliberately a bit wider than the
# active CEA cells in Tier 1 so bilinear interpolation has support).
# Tier-1 active extent: lat 28.09..43.70, lon -31.33..4.60.
# We narrow to the canonical Iberian box but extend a half-degree margin.
IBERIA_AREA = {
    "north": 44.0,
    "south": 35.0,
    "west": -10.0,
    "east": 4.0,
}

# Decade slices.
HORIZONS = {
    "2046_2055": ("20460101", "20551231"),
    "2076_2085": ("20760101", "20851231"),
}

# Minimum sane file size (in bytes) for the cached extract to be
# considered usable. Set conservatively; will skip refetch if the cached
# file is at least this big.
MIN_BYTES = 1_000_000  # 1 MB

# %% [markdown]
# ## DestinE Climate DT request template
#
# The keys below are the typical shape of a polytope request against
# the DestinE Climate DT catalogue. **Treat the param identifiers and
# stream/resolution choices as placeholders to verify on the DestinE
# platform.**

# %%
print(
    "CHECK: verify DestinE catalogue keys before first run. The keys "
    "below are placeholders. See: https://destine.ecmwf.int/ for the "
    "current Climate DT catalogue."
)


def _build_request(start_date: str, end_date: str) -> dict:
    """Polytope request body. CHECK every key on first DestinE run."""
    return {
        "class": "d1",                    # CHECK: DestinE Climate DT class
        "dataset": "climate-dt",          # CHECK: dataset identifier
        "activity": "ScenarioMIP",        # CHECK: activity (CMIP-aligned)
        "experiment": "SSP3-7.0",         # CHECK: experiment label
        "model": "ICON",                  # CHECK: model — ICON or IFS-NEMO
        "resolution": "high",             # CHECK: high vs standard
        "type": "fc",                     # CHECK: forecast type
        "stream": "clte",                 # CHECK: climate (clte) vs other
        "param": "165/166/167",           # CHECK: tmax / tmin / total precip
        "date": f"{start_date}/to/{end_date}",
        "time": "0000",
        "area": [
            IBERIA_AREA["north"],
            IBERIA_AREA["west"],
            IBERIA_AREA["south"],
            IBERIA_AREA["east"],
        ],
    }


# %% [markdown]
# ## Fetch each horizon

# %%
import earthkit.data as ekd  # noqa: E402  (deferred to keep guard cheap)

for horizon_name, (start, end) in HORIZONS.items():
    out_path = DATA_DIR / f"destine_iberia_{horizon_name}.nc"
    if out_path.exists() and out_path.stat().st_size > MIN_BYTES:
        print(f"[cached] {out_path.name}  ({out_path.stat().st_size:,} bytes)")
        continue

    request = _build_request(start, end)
    print(f"\n[fetch] {horizon_name}: {start} .. {end}")
    print(f"  request: {request}")

    # CHECK: the second positional arg is the polytope **collection**
    # (server namespace), distinct from the in-request "dataset" key.
    # For Destination Earth Climate DT the canonical collection name
    # at the time of writing is "destination-earth"; alternatives
    # encountered in the wild include "destination-earth-climate-dt"
    # and "destination-earth-data-lake". If you hit
    # `polytope.api.exceptions.UnknownCollection`, swap the value.
    POLYTOPE_COLLECTION = "destination-earth"
    print(f"  polytope collection: {POLYTOPE_COLLECTION}")

    ds = ekd.from_source("polytope", POLYTOPE_COLLECTION, request)

    # earthkit-data returns a fieldlist; convert to xarray for NetCDF I/O.
    # CHECK: the to_xarray() accessor may need keyword args (e.g.
    # ``time_dim_mode="forecast"``) depending on the request shape.
    xa = ds.to_xarray()

    print(f"  variables: {list(xa.data_vars)}")
    print(f"  coords: {list(xa.coords)}")

    xa.to_netcdf(out_path)
    print(f"  wrote {out_path}  ({out_path.stat().st_size:,} bytes)")

print("\nAll horizons present.")
for horizon_name in HORIZONS:
    p = DATA_DIR / f"destine_iberia_{horizon_name}.nc"
    print(f"  {p}  exists={p.exists()}  "
          f"size={p.stat().st_size if p.exists() else 0:,}")
