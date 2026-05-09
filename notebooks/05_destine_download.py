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
#   * **Near-term**: 2020–2029 (first decade of SSP3-7.0; analogue of
#     Soroye's "recent" 2000–2014 window)
#   * **Mid-term**: 2030–2039 (second decade — "by 2040" projection)
#
# **Why not mid-/end-of-century?** The DestinE Climate DT Phase 1
# archive is populated through 2039 inclusive (verified 2026-05-09).
# The original design called for 2046–2055 + 2076–2085 horizons; those
# are deferred to a follow-up Outcome when the archive extends past
# 2050. Document this in `nanopubs/drafts/05_outcome.md` § Limitations.
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
    "2020_2029": ("20200101", "20291231"),
    "2030_2039": ("20300101", "20391231"),
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
    # NOTE: no `area` key. DestinE IFS-NEMO output is HEALPix nside=128
    # NESTED, and MARS errors with "Representation::croppedRepresentation()
    # not implemented for HEALPixNested" if `area` is supplied. We fetch
    # globally (~786 KB per timestep, lightweight on the wire) and subset
    # to Iberian HEALPix pixels in-process before writing the NetCDF.
    return {
        "class": "d1",
        "dataset": "climate-dt",
        "activity": "ScenarioMIP",
        "experiment": "SSP3-7.0",
        "expver": "0001",
        "generation": "1",
        "realization": "1",
        "model": "IFS-NEMO",
        "resolution": "standard",
        "type": "fc",
        "stream": "clte",
        "levtype": "sfc",
        "param": param,
        "date": f"{start_date}/to/{end_date}",
        "time": time,
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

# %% [markdown]
# ## Iberia HEALPix mask (input-side only)
#
# DestinE IFS-NEMO output is on HEALPix nside=128 NESTED (~196,608
# global cells). We use **`healpix-geo`** (not `healpy`) here because
# `healpy` is non-geo-aware (cosmology-first; lon convention 0–360, no
# CRS handling, no datum) and accumulates small biases that are
# unacceptable for biodiversity-precision climate-impact work over
# Iberia (`DOMAIN.md` § Biodiversity is high-precision climate-impact
# science).
#
# **HEALPix is the INPUT format only.** This notebook subsets the
# global HEALPix output to Iberian cells and writes the subset NetCDF.
# The actual *analytical grid* is the 100 km CEA grid the upstream
# weatherxbio v0.2.0 GLMM was fitted on — the HEALPix→CEA aggregation
# happens at the import boundary in `06_destine_clean.py`. The GLMM
# coefficients are applied on CEA cells (same grid they were fitted
# on), preserving consistency with the validated upstream `sc_TEI_delta`.
#
# Per `DOMAIN.md`: HEALPix indexing is **always NESTED**.

# %%
import numpy as np  # noqa: E402
# CHECK: healpix-geo API. The exact symbol/function names may differ
# in your installed version. The minimal contract we need is:
#   - npix lookup for nside=128 NESTED  (==> 196,608 cells)
#   - geo-aware cell-center lon/lat lookup, lon in [-180, 180]
# If the installed package surfaces these under different names, swap
# the imports below; the rest of the notebook is API-agnostic.
import healpix_geo  # noqa: E402

DESTINE_NSIDE = 128
DESTINE_NPIX = 12 * DESTINE_NSIDE * DESTINE_NSIDE   # 196,608 (HEALPix invariant)

# CHECK: `healpix_geo.nested.cell_centroids(nside, ipix)` is the most
# likely API; alternatives include `healpix_geo.HealpixGrid(...).cell_centers`
# or `healpix_geo.lonlat(nside, ipix, scheme='nested')`. Adjust to match
# your version on the DestinE platform.
_all_pix = np.arange(DESTINE_NPIX)
try:
    _lons, _lats = healpix_geo.nested.cell_centroids(DESTINE_NSIDE, _all_pix)
except AttributeError:
    # Fallback shape: a HealpixGrid-class API
    _grid = healpix_geo.HealpixGrid(nside=DESTINE_NSIDE, scheme="nested")
    _lons, _lats = _grid.cell_centers(_all_pix)

# healpix-geo is geo-aware → lon in [-180, 180], no Greenwich-wrap kludge
_iberia_mask = (
    (_lons >= IBERIA_AREA["west"])  & (_lons <= IBERIA_AREA["east"])
    & (_lats >= IBERIA_AREA["south"]) & (_lats <= IBERIA_AREA["north"])
)
IBERIA_PIX = _all_pix[_iberia_mask]

print(f"Iberia HEALPix nside={DESTINE_NSIDE} NESTED: "
      f"{len(IBERIA_PIX):,} cells / {DESTINE_NPIX:,} global "
      f"({100 * len(IBERIA_PIX) / DESTINE_NPIX:.2f}%)")
print("HEALPix is INPUT format only — analytical grid is CEA "
      "(aggregation happens in 06_destine_clean.py).")

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
        print(f"  global sizes: {dict(xa.sizes)}")

        # Subset to Iberian HEALPix cells before writing — DestinE
        # redistribution licence allows aggregated derived statistics
        # but not the global raw extract.
        # CHECK: the spatial dim name on DestinE is usually "values"
        # (HEALPix payload). If a different name is exposed, swap below.
        spatial_dim = "values" if "values" in xa.sizes else (
            "cell" if "cell" in xa.sizes else None
        )
        if spatial_dim is None:
            raise SystemExit(
                "Could not find HEALPix spatial dim in xarray output. "
                f"Available dims: {dict(xa.sizes)}"
            )
        if xa.sizes[spatial_dim] != DESTINE_NPIX:
            print(f"  WARNING: spatial dim '{spatial_dim}' has "
                  f"{xa.sizes[spatial_dim]:,} cells, expected "
                  f"{DESTINE_NPIX:,} for nside={DESTINE_NSIDE}. "
                  "Iberia subset may be wrong — verify HEALPix nside.")

        xa_iberia = xa.isel({spatial_dim: IBERIA_PIX})
        # Stash Iberia HEALPix indices as a coord so 06 can rebuild
        # cell positions without reapplying the bbox mask.
        xa_iberia = xa_iberia.assign_coords(
            iberia_pix=(spatial_dim, IBERIA_PIX)
        )
        print(f"  Iberia subset sizes: {dict(xa_iberia.sizes)}")

        xa_iberia.to_netcdf(out_path)
        print(f"  wrote {out_path}  ({out_path.stat().st_size:,} bytes)")

print("\nAll (horizon × variable) extracts present.")
for horizon_name in HORIZONS:
    for spec in VARIABLE_SPECS:
        p = DATA_DIR / f"destine_iberia_{horizon_name}_{spec['label']}.nc"
        print(f"  {p.name}  exists={p.exists()}  "
              f"size={p.stat().st_size if p.exists() else 0:,}")
