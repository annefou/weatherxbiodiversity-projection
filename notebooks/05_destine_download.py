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
# # 05 — DestinE Climate DT GRIB retrieve (DestinE platform only)
#
# **Minimal-on-DestinE design**: this notebook does only the polytope
# retrieval. It writes the raw global GRIB files to
# `data/destine/raw/destine_<horizon>_<var>.grib`. **No GRIB decoding,
# no Iberia subset, no xarray** runs here — those steps need
# `eccodes` + a HEALPix-NESTED-aware decode path that the DestinE
# platform's eccodes (RING-only Geoiterator) cannot provide.
#
# After this notebook completes:
#
# 1. Tar the four GRIBs and download to local Mac (instructions
#    printed at the end of this notebook).
# 2. On local Mac, run `notebooks/06_destine_clean.py` which decodes
#    the GRIBs (eccodes Python API, no Geoiterator), subsets to the
#    pre-computed Iberian HEALPix-NESTED nside=128 cells, aggregates
#    via NESTED parent operation to nside=64 (Tier-1 analytical grid),
#    and computes future TEI/PEI per species per cell.
#
# ## SSP3-7.0 horizons
#
# The DestinE Climate DT Phase 1 archive is currently populated through
# 2039 inclusive (verified 2026-05-09). Horizons:
#
#   * **Near-term**: 2020–2029
#   * **Mid-term**:  2030–2039
#
# Mid-/end-of-century horizons (2046–2055 + 2076–2085) are deferred to
# a follow-up Outcome when the archive extends past 2050.
#
# ## Variables
#
# Two requests per horizon (instantaneous t2m + accumulated tp can't be
# combined in one polytope request):
#
#   * **`param=167`** (2m temperature, K, instantaneous, 4×/day at
#     00/06/12/18 UTC) — daily max/min derived in 06.
#   * **`param=228`** (total precipitation, m, accumulated, 1×/day at
#     0000 UTC) — 24-hour daily totals.
#
# Total = 4 polytope retrievals (2 horizons × 2 variables).

# %%
import tempfile
from pathlib import Path

from polytope.api import Client

# %%
ROOT = Path("..").resolve()
GRIB_DIR = ROOT / "data" / "destine" / "raw"
GRIB_DIR.mkdir(parents=True, exist_ok=True)

POLYTOPE_COLLECTION = "destination-earth"
POLYTOPE_ADDRESS = "https://polytope.lumi.apps.dte.destination-earth.eu"

# Decade slices (verified populated in the Climate DT archive).
HORIZONS = {
    "2020_2029": ("20200101", "20291231"),
    "2030_2039": ("20300101", "20391231"),
}

VARIABLE_SPECS = [
    {"label": "t2m", "param": "167",
     "time": "0000/0600/1200/1800", "encoding": "instantaneous"},
    {"label": "tp",  "param": "228",
     "time": "0000",                "encoding": "accumulated"},
]

# 1 MB cache threshold (any bigger = treat as a usable retrieval).
MIN_BYTES = 1_000_000


# %%
def build_request(start_date: str, end_date: str, *,
                  param: str, time: str) -> dict:
    """Polytope request body for DestinE Climate DT SSP3-7.0
    IFS-NEMO at 'standard' resolution (HEALPix nside=128 NESTED).
    Verified-working keys per the DestinE platform documentation
    (2026-05-09): expver='0001', generation='1', realization='1'.

    No `area` key — DestinE IFS-NEMO is HEALPix-archived and MARS
    can't crop HEALPix data with a lat/lon bbox (raises
    "Representation::croppedRepresentation() not implemented for
    HEALPixNested"). We fetch globally and subset Iberia in 06.
    """
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


# %% [markdown]
# ## Fetch each (horizon × variable)

# %%
client = Client(address=POLYTOPE_ADDRESS)

print(f"polytope collection: {POLYTOPE_COLLECTION}")
print(f"polytope address:    {POLYTOPE_ADDRESS}\n")

for horizon_name, (start, end) in HORIZONS.items():
    for spec in VARIABLE_SPECS:
        out_grib = GRIB_DIR / (
            f"destine_{horizon_name}_{spec['label']}.grib"
        )
        if out_grib.exists() and out_grib.stat().st_size > MIN_BYTES:
            print(f"[cached] {out_grib.name}  "
                  f"({out_grib.stat().st_size:,} bytes)")
            continue

        request = build_request(
            start, end,
            param=spec["param"],
            time=spec["time"],
        )
        print(f"\n[fetch] {horizon_name} / {spec['label']} "
              f"({spec['encoding']}): {start} .. {end}")
        print(f"  request: {request}")

        client.retrieve(
            POLYTOPE_COLLECTION, request,
            output_file=str(out_grib),
            asynchronous=False,
        )
        print(f"  saved → {out_grib} "
              f"({out_grib.stat().st_size:,} bytes)")

# %% [markdown]
# ## Final inventory + transfer instructions

# %%
print("\nGRIB inventory:")
total = 0
for horizon_name in HORIZONS:
    for spec in VARIABLE_SPECS:
        p = GRIB_DIR / f"destine_{horizon_name}_{spec['label']}.grib"
        sz = p.stat().st_size if p.exists() else 0
        total += sz
        print(f"  {p.name}  exists={p.exists()}  "
              f"size={sz:>15,} bytes ({sz / 1e9:.2f} GB)")
print(f"  Total: {total:,} bytes ({total / 1e9:.2f} GB)")

print(
    "\nNext steps (transfer GRIBs to local Mac, then run 06 there):\n"
    "  1. tar czf /tmp/destine_gribs.tar.gz -C ../data/destine raw/\n"
    "  2. Download /tmp/destine_gribs.tar.gz via JupyterLab file browser.\n"
    "  3. On local Mac, in repo root:\n"
    "       mkdir -p data/destine && \\\n"
    "       tar xzf ~/Downloads/destine_gribs.tar.gz -C data/destine/\n"
    "  4. Run snakemake locally — 06 (clean) will pick up the GRIBs."
)
