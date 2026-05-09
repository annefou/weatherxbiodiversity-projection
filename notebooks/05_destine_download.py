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
# ## DestinE Climate DT request — known constraints
#
# Verified against the official catalogue
# (`DestinE ClimateDT Parameters - DGOV` Confluence page, retrieved
# 2026-05). Two findings shape this request:
#
# 1. **DestinE Climate DT has NO daily max/min temperature encodings.**
#    Sections 4 and 5 of the catalogue explicitly state "There are no
#    maximum/minimum encodings for this dataset." We therefore fetch
#    hourly 2m temperature (`param=167`, instantaneous, levtype=sfc)
#    and derive daily max/min in `06_destine_clean.py`. To keep the
#    extract size manageable we sample 4 times/day
#    (00, 06, 12, 18 UTC) — that captures the diurnal range well
#    enough for a Soroye-style monthly-statistics analysis.
# 2. **Total precipitation (`param=228`) lives in section 3
#    (accumulated encodings), distinct from instantaneous fields.**
#    Mixing instantaneous (167) and accumulated (228) in one
#    polytope request is brittle on MARS — we issue them as TWO
#    separate requests below.
#
# All values below are catalogue-verified (model=ICON, levtype=sfc,
# stream=clte, type=fc) for the SSP3-7.0 ScenarioMIP run.

# %%
def _build_request(start_date: str, end_date: str, *,
                   param: str, time: str, encoding: str) -> dict:
    """Polytope request body, aligned with the verified-working
    DestinE Climate DT request template provided in the platform
    documentation. ``encoding`` is informational only — the section
    of the catalogue (instantaneous / accumulated) we're querying.

    Required keys for SSP3-7.0 retrieval (verified 2026-05-09):
      expver='0001', generation='1', realization='1',
      model='IFS-NEMO', resolution='standard'.
    Without these the request expands to the right cardinality but
    matches no archived data (MARS returns "0 messages retrieved").
    """
    return {
        "class": "d1",                          # DestinE
        "dataset": "climate-dt",                # Climate DT
        "activity": "ScenarioMIP",
        "experiment": "SSP3-7.0",
        "expver": "0001",                       # experiment version (required)
        "generation": "1",                      # model generation (required)
        "realization": "1",                     # ensemble member (required)
        "model": "IFS-NEMO",                    # IFS-NEMO has SSP3-7.0; ICON did not match
        "resolution": "standard",               # 'standard' is the archived resolution
        "type": "fc",
        "stream": "clte",                       # Climate experimental
        "levtype": "sfc",                       # surface fields
        "param": param,                         # 167 (2t) or 228 (tp)
        "date": f"{start_date}/to/{end_date}",
        "time": time,                           # 4-times/day for 167; '0000' for 228
        "area": [
            IBERIA_AREA["north"],
            IBERIA_AREA["west"],
            IBERIA_AREA["south"],
            IBERIA_AREA["east"],
        ],
    }


# Per-variable request specs. Each horizon is fetched twice — once for
# 2m temperature (instantaneous, 4 samples/day → derive daily max/min
# in 06_destine_clean.py), once for total precipitation
# (accumulated, daily 0000Z value gives 24-hour accumulation).
VARIABLE_SPECS = [
    {"param": "167", "time": "0000/0600/1200/1800",
     "encoding": "instantaneous", "label": "t2m"},
    {"param": "228", "time": "0000",
     "encoding": "accumulated",   "label": "tp"},
]


# %% [markdown]
# ## Fetch each horizon

# %%
import earthkit.data as ekd  # noqa: E402  (deferred to keep guard cheap)

# Polytope server (LUMI-hosted DestinE) and collection ("destination-earth").
# DestinE-issued tokens are NOT accepted by ECMWF's general
# https://polytope.ecmwf.int — must use the LUMI URL below.
POLYTOPE_COLLECTION = "destination-earth"
POLYTOPE_ADDRESS = "https://polytope.lumi.apps.dte.destination-earth.eu"
# Alternatives: "https://polytope.destination-earth.eu"
#               "https://polytope-climate-dt.destination-earth.eu"

print(f"polytope collection: {POLYTOPE_COLLECTION}")
print(f"polytope address:    {POLYTOPE_ADDRESS}\n")

for horizon_name, (start, end) in HORIZONS.items():
    for spec in VARIABLE_SPECS:
        out_path = DATA_DIR / (
            f"destine_iberia_{horizon_name}_{spec['label']}.nc"
        )
        if out_path.exists() and out_path.stat().st_size > MIN_BYTES:
            print(f"[cached] {out_path.name}  "
                  f"({out_path.stat().st_size:,} bytes)")
            continue

        request = _build_request(
            start, end,
            param=spec["param"],
            time=spec["time"],
            encoding=spec["encoding"],
        )
        print(f"\n[fetch] {horizon_name} / {spec['label']} "
              f"({spec['encoding']}): {start} .. {end}")
        print(f"  request: {request}")

        ds = ekd.from_source(
            "polytope",
            POLYTOPE_COLLECTION,
            request,
            address=POLYTOPE_ADDRESS,
        )

        # earthkit-data returns a fieldlist; convert to xarray for NetCDF I/O.
        # CHECK: the to_xarray() accessor may need keyword args (e.g.
        # ``time_dim_mode="forecast"``) depending on the request shape.
        xa = ds.to_xarray()

        print(f"  variables: {list(xa.data_vars)}")
        print(f"  coords: {list(xa.coords)}")

        xa.to_netcdf(out_path)
        print(f"  wrote {out_path}  ({out_path.stat().st_size:,} bytes)")

print("\nAll (horizon × variable) extracts present.")
for horizon_name in HORIZONS:
    for spec in VARIABLE_SPECS:
        p = DATA_DIR / f"destine_iberia_{horizon_name}_{spec['label']}.nc"
        print(f"  {p.name}  exists={p.exists()}  "
              f"size={p.stat().st_size if p.exists() else 0:,}")
