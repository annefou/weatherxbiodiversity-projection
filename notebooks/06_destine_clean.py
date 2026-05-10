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
# # 06 — DestinE clean: future TEI / PEI on the HEALPix-NESTED nside=64 substrate
#
# **Tier-2 Phase D entry point** — runs *off-DestinE* (i.e. on the local
# Mac after the four raw GRIBs have been fetched + transferred via
# `notebooks/05_destine_download.py`). The historical `_tier2_guard`
# import is deliberately omitted: 06/07/08 only need the GRIBs to be
# present on disk; if they are not, Snakemake's input-not-found error
# is a clearer signal than a silent skip.
#
# Pipeline:
#
#   1. Decode each raw GRIB **per message** via the eccodes Python API
#      (no cfgrib, no Geoiterator — the DestinE GRIBs are HEALPix
#      NESTED Nside=128 and eccodes' Geoiterator is RING-only).
#   2. Subset each message to the 440 pre-computed Iberian nside=128
#      NESTED cells immediately at decode time, keeping memory bounded
#      to (n_msgs × 440 × 8 bytes).
#   3. Aggregate timestepwise to monthly per nside=128 cell:
#        * `t2m` (4 instantaneous samples/day at 00/06/12/18 UTC) →
#          per-day max + per-day min across the 4 samples → monthly
#          max-of-daily-max and monthly min-of-daily-min per cell. This
#          is the standard ECMWF "daily Tmax/Tmin from 6-hourly
#          snapshots" approximation; we also derive a daily-mean from
#          the 4 samples and a monthly mean of daily means for the CRU
#          'tmp' analogue.
#        * `tp` (1 accumulated value/day at step=1) → monthly sum of
#          daily totals (in metres).
#   4. Aggregate nside=128 → nside=64 via the NESTED parent operation
#      (4-fold equal-area mean): `parent = child >> 2`. Built from the
#      precomputed (110 × 4) lookup table; arithmetic only — no
#      resampling, no healpix_geo dependency in this step.
#   5. Compute future-decade mean climate per nside=64 cell:
#        * `meanT_future[c]` = monthly-tmean averaged across the decade
#          (matches Tier-1 CRU `tmp_annual` mean-of-monthly-means).
#        * `meanP_future[c]` = annual total precipitation averaged
#          across the decade (mm/yr; matches Tier-1 CRU `pre_annual`
#          which is sum of monthly totals → period mean of annual sum).
#   6. Compute future TEI / PEI per species per nside=64 cell using
#      historical species niche limits (`T_min_spp`, `T_max_spp`,
#      `P_min_spp`, `P_max_spp`) loaded from
#      `healpix_port/outputs_iberia/climate_tei_pei_healpix.npz`. The
#      historical TEI_bs / PEI_bs are also carried verbatim so 07 has
#      both `*_bs` and `*_delta` predictors with no recomputation.
#
# Outputs (gitignored — see `.gitignore` Tier-2 section):
#
#   * `healpix_port/outputs_iberia/climate_tei_pei_future_2020_2029_healpix.nc`
#   * `healpix_port/outputs_iberia/climate_tei_pei_future_2030_2039_healpix.nc`

# %%
from datetime import date
from pathlib import Path

import eccodes
import numpy as np
import pandas as pd
import xarray as xr

# Make the healpix_port package importable when run with cwd=notebooks/
import sys as _sys
_sys.path.insert(0, str(Path("..").resolve()))
from healpix_port._dggs_metadata import PROJECT_DGGS_ATTRS  # noqa: E402

# %%
ROOT = Path("..").resolve()
GRIB_DIR = ROOT / "data" / "destine" / "raw"
HPORT = ROOT / "healpix_port"
OUT_DIR = HPORT / "outputs_iberia"
PRECOMP = ROOT / "data" / "precomputed"

HORIZONS = ["2020_2029", "2030_2039"]

# DestinE Climate DT clte-stream tp is HOURLY accumulation; the
# 05_destine_download.py request fetched only `time=0000`, leaving
# 1 hour per day. Multiply by 24 to upscale to a climatological
# daily-total estimate. See docstring of `aggregate_tp` for details.
TP_HOURLY_TO_DAILY_FACTOR = 24

# %% [markdown]
# ## Iberian HEALPix-NESTED indexing
#
# `iberia_pix_nside128_nested.npy` holds the 440 nside=128 cells we
# subset to during decode. `iberia_pix_nside64_nested.npy` holds the
# 110 parent cells that match the Tier-1 analytical grid. By
# construction (see `scripts/precompute_iberia_pix.py`) every nside=128
# child has its parent in the nside=64 list — the parent operation is
# the bit-shift `child >> 2`.

# %%
IBERIA_PIX_128 = np.load(PRECOMP / "iberia_pix_nside128_nested.npy").astype(np.int64)
IBERIA_PIX_64 = np.load(PRECOMP / "iberia_pix_nside64_nested.npy").astype(np.int64)
N_128 = len(IBERIA_PIX_128)
N_64 = len(IBERIA_PIX_64)
print(f"Iberia substrate: {N_64} nside=64 parents × 4 = {N_128} nside=128 children")

# Build the (N_64, 4) child-index lookup: for each parent (row), the
# four positions in IBERIA_PIX_128 whose parent matches. Used for the
# equal-area parent-mean aggregation.
parents_of_children = (IBERIA_PIX_128 >> 2)
child_pos_for_parent = np.full((N_64, 4), -1, dtype=np.int64)
counter = np.zeros(N_64, dtype=np.int64)
parent_to_row = {int(p): i for i, p in enumerate(IBERIA_PIX_64)}
for j, parent in enumerate(parents_of_children):
    row = parent_to_row[int(parent)]
    child_pos_for_parent[row, counter[row]] = j
    counter[row] += 1
assert (counter == 4).all(), "every parent must have exactly 4 children"
print("Built (110, 4) child-index lookup for parent-mean aggregation.")

# Inverse lookup: for each position j in IBERIA_PIX_128, the row index
# of its parent in IBERIA_PIX_64. Used by the Option-B nside=128
# projection path which inherits historical baselines from the
# nside=64 GLMM training data per parent.
parent_row_for_child = np.array(
    [parent_to_row[int(p)] for p in parents_of_children], dtype=np.int64
)
print(f"Built ({N_128},) parent-row lookup for nside=128 child-to-parent inheritance.")

# %% [markdown]
# ## Historical species niche limits (Tier-1 HEALPix), held fixed

# %%
hist = xr.open_dataset(OUT_DIR / "climate_tei_pei_healpix.nc")
species = [str(s) for s in hist["species"].values]
n_spp = len(species)
T_min_spp = hist["T_min_spp"].values
T_max_spp = hist["T_max_spp"].values
P_min_spp = hist["P_min_spp"].values
P_max_spp = hist["P_max_spp"].values
TEI_bs_hist = hist["tei_bs"].values        # (n_spp, 110)
PEI_bs_hist = hist["pei_bs"].values
avgtemp_bs_hist = hist["avgtemp_bs"].values   # (110,)
avgprecip_bs_hist = hist["avgprecip_bs"].values
iberia_cells_hp = hist["cell_ids"].values.astype(np.uint64)
print(f"Loaded historical niche limits: {n_spp} species, {len(iberia_cells_hp)} cells")

# Sanity: iberia_cells_hp must equal IBERIA_PIX_64 (Tier-1 was fitted
# on these). If not, the parent-aggregation alignment is wrong.
assert np.array_equal(iberia_cells_hp.astype(np.int64), IBERIA_PIX_64), (
    "Tier-1 HEALPix substrate cells do not match the precomputed nside=64 list."
)


# %% [markdown]
# ## Per-message GRIB decoder
#
# Streams every message in the GRIB and yields
# `(timestamp, values_iberia_128)`. We bypass cfgrib entirely and
# never touch latitudes/longitudes or the Geoiterator — the
# orderingConvention='nested' guarantee + the precomputed Iberian
# nside=128 index means we can just do `values[IBERIA_PIX_128]` to
# subset.

# %%
def stream_grib(grib_path: Path):
    """Yield (timestamp, values_iberia_128, dataDate_int, step_int) for every message."""
    with open(grib_path, "rb") as f:
        msg_count = 0
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                # Defensive sanity-checks on the first message only.
                if msg_count == 0:
                    if eccodes.codes_get(gid, "Nside") != 128:
                        raise RuntimeError(
                            f"Expected Nside=128, got {eccodes.codes_get(gid, 'Nside')}"
                        )
                    if eccodes.codes_get(gid, "orderingConvention") != "nested":
                        raise RuntimeError(
                            "Expected orderingConvention='nested' "
                            f"got {eccodes.codes_get(gid, 'orderingConvention')!r}"
                        )

                date = eccodes.codes_get(gid, "dataDate")  # int YYYYMMDD
                time = eccodes.codes_get(gid, "dataTime")  # int e.g. 0, 600, 1200, 1800
                step = eccodes.codes_get(gid, "step")
                values = eccodes.codes_get_array(gid, "values")
                # Subset to Iberia nside=128 (440 cells).
                vals_iberia = values[IBERIA_PIX_128].astype(np.float64, copy=True)
            finally:
                eccodes.codes_release(gid)

            hour = time // 100
            minute = time % 100
            ts = pd.Timestamp(
                year=date // 10000,
                month=(date // 100) % 100,
                day=date % 100,
                hour=int(hour),
                minute=int(minute),
            )
            yield ts, vals_iberia, int(date), int(step)
            msg_count += 1


# %% [markdown]
# ## Per-horizon driver

# %%
def aggregate_t2m(grib_path: Path):
    """Decode the 4×/day instantaneous t2m GRIB.

    Returns three (n_months, N_128) arrays of monthly Tmax, Tmin, Tmean
    per Iberian nside=128 cell. Tmax = monthly max-of-daily-max,
    Tmin = monthly min-of-daily-min, Tmean = monthly mean-of-daily-mean.
    """
    print(f"  decoding {grib_path.name} ...")
    # Keep per-day max / min / running mean per cell, indexed by date.
    # Use float32 dicts to bound memory (each value: 440 cells × 8 B).
    daily_max: dict[pd.Timestamp, np.ndarray] = {}
    daily_min: dict[pd.Timestamp, np.ndarray] = {}
    daily_sum: dict[pd.Timestamp, np.ndarray] = {}
    daily_n: dict[pd.Timestamp, int] = {}

    n_msgs = 0
    for ts, vals, date, step in stream_grib(grib_path):
        # Validity day = the calendar day of the message (step=0 instant).
        day = ts.normalize()
        if day in daily_max:
            np.maximum(daily_max[day], vals, out=daily_max[day])
            np.minimum(daily_min[day], vals, out=daily_min[day])
            daily_sum[day] += vals
            daily_n[day] += 1
        else:
            daily_max[day] = vals.copy()
            daily_min[day] = vals.copy()
            daily_sum[day] = vals.copy()
            daily_n[day] = 1
        n_msgs += 1
        if n_msgs % 2000 == 0:
            print(f"    {n_msgs:>6} messages, {len(daily_max):>5} days so far")
    print(f"    decoded {n_msgs} messages → {len(daily_max)} days")

    # Stack daily arrays into (n_days, N_128).
    days = sorted(daily_max.keys())
    dmax = np.stack([daily_max[d] for d in days], axis=0)        # (n_days, 440)
    dmin = np.stack([daily_min[d] for d in days], axis=0)
    dsum = np.stack([daily_sum[d] for d in days], axis=0)
    dn = np.array([daily_n[d] for d in days])[:, None]           # (n_days, 1)
    dmean = dsum / dn

    # Monthly aggregation: monthly max-of-daily-max, etc. We use a
    # pandas DatetimeIndex grouped by year-month for the boundaries.
    days_idx = pd.DatetimeIndex(days)
    months = days_idx.to_period("M")
    unique_months = months.unique()
    n_months = len(unique_months)
    mmax = np.empty((n_months, N_128), dtype=np.float64)
    mmin = np.empty((n_months, N_128), dtype=np.float64)
    mmean = np.empty((n_months, N_128), dtype=np.float64)
    for i, m in enumerate(unique_months):
        mask = months == m
        mmax[i] = dmax[mask].max(axis=0)
        mmin[i] = dmin[mask].min(axis=0)
        mmean[i] = dmean[mask].mean(axis=0)
    print(f"    monthly arrays: {n_months} months × {N_128} cells")
    return unique_months, mmax, mmin, mmean


def aggregate_tp(grib_path: Path):
    """Decode the 1×/day accumulated tp GRIB.

    Returns (n_months, N_128) of monthly precipitation totals (in
    metres) per Iberian nside=128 cell.

    NOTE — DestinE Climate DT clte-stream tp is **hourly** accumulation
    (verified by `lengthOfTimeRange=1, indicatorOfUnitForTimeRange=1
    [hours], stepRange='0-1'`). The 05_destine_download.py request
    asked for `time=0000` only, so we have 1 message per day = 1 hour
    of accumulation centred on 0000 UTC, NOT a full daily total.
    Without correction the decadal Iberia mean is ~22 mm/yr (vs the
    CRU TS historical baseline of ~492 mm/yr) — a clear factor-of-24
    underestimate. We apply `TP_HOURLY_TO_DAILY_FACTOR = 24` to
    upscale to a climatological daily-total estimate; this assumes
    the 0000 UTC hour's precipitation rate is representative of the
    full 24-hour mean rate, which is a known approximation for
    decade-scale period means but would be invalid for sub-monthly
    statistics. Future work: re-retrieve with `time=0000/0100/.../2300`
    and sum the 24 hourly accumulations per day.
    """
    print(f"  decoding {grib_path.name} ...")
    daily_total: dict[pd.Timestamp, np.ndarray] = {}
    n_msgs = 0
    for ts, vals, date, step in stream_grib(grib_path):
        day = ts.normalize()
        # By the DestinE schema each tp message = 1 day total at step=1.
        # Defensive: if multiple messages share a day, sum (the cumulative
        # convention can vary across runs).
        if day in daily_total:
            daily_total[day] += vals
        else:
            daily_total[day] = vals.copy()
        n_msgs += 1
        if n_msgs % 1000 == 0:
            print(f"    {n_msgs:>6} messages, {len(daily_total):>5} days so far")
    print(f"    decoded {n_msgs} messages → {len(daily_total)} days")

    days = sorted(daily_total.keys())
    dtot = np.stack([daily_total[d] for d in days], axis=0)        # (n_days, 440)
    days_idx = pd.DatetimeIndex(days)
    months = days_idx.to_period("M")
    unique_months = months.unique()
    n_months = len(unique_months)
    msum = np.empty((n_months, N_128), dtype=np.float64)
    for i, m in enumerate(unique_months):
        mask = months == m
        msum[i] = dtot[mask].sum(axis=0)
    print(f"    monthly precip arrays: {n_months} months × {N_128} cells")
    return unique_months, msum


def parent_mean_128_to_64(arr_128: np.ndarray) -> np.ndarray:
    """Equal-area mean of nside=128 → nside=64 via the NESTED parent op.

    Input: arr_128 of shape (..., N_128).
    Output: arr_64 of shape (..., N_64).
    Uses the precomputed (N_64, 4) child-position lookup.
    """
    children = arr_128[..., child_pos_for_parent]   # (..., N_64, 4)
    return children.mean(axis=-1)


# %% [markdown]
# ## Run aggregation per horizon

# %%
for horizon in HORIZONS:
    print(f"\n=== Horizon {horizon} ===")
    t2m_path = GRIB_DIR / f"destine_{horizon}_t2m.grib"
    tp_path = GRIB_DIR / f"destine_{horizon}_tp.grib"
    if not (t2m_path.exists() and tp_path.exists()):
        print(f"  [skip] missing GRIBs for {horizon}")
        continue

    months_t, mmax128, mmin128, mmean128 = aggregate_t2m(t2m_path)
    months_p, msum128 = aggregate_tp(tp_path)

    # Convert temperatures K → degC (Tier-1 CRU is in degC).
    mmax128_c = mmax128 - 273.15
    mmin128_c = mmin128 - 273.15
    mmean128_c = mmean128 - 273.15

    # Aggregate to nside=64 via NESTED parent equal-area mean.
    mmean64 = parent_mean_128_to_64(mmean128_c)        # (n_months, N_64)
    mmax64 = parent_mean_128_to_64(mmax128_c)
    mmin64 = parent_mean_128_to_64(mmin128_c)
    msum64 = parent_mean_128_to_64(msum128)            # (n_months, N_64) (m/month)

    # ---- Tier-1-matched per-cell summary statistics ----
    # Tier-1 (healpix_port/04_climate_tei_pei_healpix.py):
    #   tmp_annual = monthly tmp resampled annually with mean()
    #     -> meanT_bs = nanmean of annual values across years
    #   pre_annual = monthly pre resampled annually with sum()
    #     -> meanP_bs = nanmean of annual sums across years (mm/yr)
    # Future contract: same definitions, 10-year decade as the
    # "period". meanT_future = mean over months over years = mean of
    # all monthly tmean. meanP_future = mean of annual totals.

    meanT_future_64 = mmean64.mean(axis=0)            # (N_64,)  degC

    # Annual precip totals: sum monthly within each year, then mean
    # across years. Build year index from `months_p`.
    years_p = np.array([m.year for m in months_p])
    unique_years = np.unique(years_p)
    annual_p_64 = np.empty((len(unique_years), N_64))
    for i, y in enumerate(unique_years):
        sel = years_p == y
        annual_p_64[i] = msum64[sel].sum(axis=0)
    # Convert to mm (Tier-1 CRU pre is in mm; DestinE tp is metres)
    # AND apply the hourly→daily upscale factor (see TP note above).
    annual_p_64_mm = annual_p_64 * 1000.0 * TP_HOURLY_TO_DAILY_FACTOR
    meanP_future_64 = annual_p_64_mm.mean(axis=0)     # (N_64,) mm/yr

    print(
        f"  meanT_future:  {meanT_future_64.min():.2f} .. "
        f"{meanT_future_64.max():.2f} degC  "
        f"(median {np.median(meanT_future_64):.2f})"
    )
    print(
        f"  meanP_future:  {meanP_future_64.min():.0f} .. "
        f"{meanP_future_64.max():.0f} mm/yr  "
        f"(median {np.median(meanP_future_64):.0f})"
    )
    print(
        f"  avgtemp_bs (Tier-1 historical) median: "
        f"{np.median(avgtemp_bs_hist):.2f} degC -- "
        f"avgtemp_delta median: "
        f"{np.median(meanT_future_64 - avgtemp_bs_hist):.2f} degC"
    )

    # Future TEI / PEI per species per nside=64 cell.
    T_range = T_max_spp - T_min_spp
    P_range = P_max_spp - P_min_spp
    with np.errstate(invalid="ignore", divide="ignore"):
        TEI_future = (
            (meanT_future_64[np.newaxis, :] - T_min_spp[:, np.newaxis])
            / T_range[:, np.newaxis]
        )
        PEI_future = (
            (meanP_future_64[np.newaxis, :] - P_min_spp[:, np.newaxis])
            / P_range[:, np.newaxis]
        )

    TEI_delta_future = (TEI_future - TEI_bs_hist).astype(np.float32)
    PEI_delta_future = (PEI_future - PEI_bs_hist).astype(np.float32)

    avgtemp_delta_future = (meanT_future_64 - avgtemp_bs_hist).astype(np.float32)
    avgprecip_delta_future = (meanP_future_64 - avgprecip_bs_hist).astype(np.float32)

    # Save as CF-compliant NetCDF — same variable layout as Tier-1
    # historical climate_tei_pei_healpix.nc so 07 can swap directly.
    out_path = OUT_DIR / f"climate_tei_pei_future_{horizon}_healpix.nc"
    decade_label = horizon.replace("_", "-")
    ds_out = xr.Dataset(
        data_vars={
            # Tier-1 historical baselines carried verbatim:
            "tei_bs": (
                ('species', 'cells'),
                TEI_bs_hist.astype(np.float32),
                {
                    "long_name": "Climatic Position Index (thermal) baseline 1901-1974",
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                    "comment": "Carried from Tier-1 historical climate_tei_pei_healpix.nc",
                },
            ),
            "pei_bs": (
                ('species', 'cells'),
                PEI_bs_hist.astype(np.float32),
                {
                    "long_name": "Climatic Position Index (precipitation) baseline 1901-1974",
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                    "comment": "Carried from Tier-1 historical climate_tei_pei_healpix.nc",
                },
            ),
            "avgtemp_bs": (
                ('cells',),
                avgtemp_bs_hist.astype(np.float32),
                {
                    "long_name": "baseline mean annual temperature",
                    "units": "degC",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "avgprecip_bs": (
                ('cells',),
                avgprecip_bs_hist.astype(np.float32),
                {
                    "long_name": "baseline mean annual total precipitation",
                    "units": "mm",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            # DestinE-derived future climate:
            "tei_future": (
                ('species', 'cells'),
                TEI_future.astype(np.float32),
                {
                    "long_name": (
                        f"Climatic Position Index (thermal) future decade {decade_label}"
                    ),
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "pei_future": (
                ('species', 'cells'),
                PEI_future.astype(np.float32),
                {
                    "long_name": (
                        f"Climatic Position Index (precipitation) future decade {decade_label}"
                    ),
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "tei_delta": (
                ('species', 'cells'),
                TEI_delta_future,
                {
                    "long_name": "Delta thermal CPI (future minus baseline)",
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "pei_delta": (
                ('species', 'cells'),
                PEI_delta_future,
                {
                    "long_name": "Delta precipitation CPI (future minus baseline)",
                    "units": "1",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "avgtemp_delta": (
                ('cells',),
                avgtemp_delta_future,
                {
                    "long_name": "delta annual mean temperature (future minus baseline)",
                    "units": "degC",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "avgprecip_delta": (
                ('cells',),
                avgprecip_delta_future,
                {
                    "long_name": "delta annual total precipitation (future minus baseline)",
                    "units": "mm",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "meanT_future": (
                ('cells',),
                meanT_future_64.astype(np.float32),
                {
                    "long_name": (
                        f"future-decade ({decade_label}) mean annual temperature"
                    ),
                    "units": "degC",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "meanP_future": (
                ('cells',),
                meanP_future_64.astype(np.float32),
                {
                    "long_name": (
                        f"future-decade ({decade_label}) mean annual total precipitation"
                    ),
                    "units": "mm",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "T_min_spp": (
                ("species",),
                T_min_spp.astype(np.float32),
                {
                    "long_name": "species-specific cold thermal limit",
                    "units": "degC",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "T_max_spp": (
                ("species",),
                T_max_spp.astype(np.float32),
                {
                    "long_name": "species-specific hot thermal limit",
                    "units": "degC",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "P_min_spp": (
                ("species",),
                P_min_spp.astype(np.float32),
                {
                    "long_name": "species-specific dry-period precipitation limit",
                    "units": "mm",
                    "_FillValue": np.float32(np.nan),
                },
            ),
            "P_max_spp": (
                ("species",),
                P_max_spp.astype(np.float32),
                {
                    "long_name": "species-specific wet-period precipitation limit",
                    "units": "mm",
                    "_FillValue": np.float32(np.nan),
                },
            ),
        },
        coords={
            "species": np.array(species, dtype=object),
            "cell_ids": ("cells", IBERIA_PIX_64.astype(np.int64)),        },
        attrs={
            "Conventions": "CF-1.10",
            "title": (
                f"Future TEI / PEI per species per cell ({decade_label}) on "
                "Iberian HEALPix nside=64 NESTED"
            ),
            "horizon": decade_label,
            "source": (
                "DestinE Climate DT SSP3-7.0 IFS-NEMO standard "
                "(HEALPix nside=128 NESTED)"
            ),
            "license_note": (
                "DestinE Climate DT redistribution restricted; aggregated "
                "derived statistics only"
            ),
            "processing": (
                "06_destine_clean.py: eccodes-direct decode -> IBERIA_PIX_128 "
                "subset -> NESTED parent (>>2) 4:1 mean -> monthly aggregation "
                "-> CPI"
            ),
            "history": (
                f"Created {date.today().isoformat()} by "
                "notebooks/06_destine_clean.py"
            ),
            **PROJECT_DGGS_ATTRS,    # DGGS Zarr Convention v1 — see healpix_port/_dggs_metadata.py
            "tp_hourly_to_daily_factor": TP_HOURLY_TO_DAILY_FACTOR,
            "n_cells": int(N_64),
        },
    )
    ds_out["cell_ids"].attrs.update({
        "long_name": "HEALPix NESTED pixel index (nside=64)",
    })
    ds_out["species"].attrs.update({"long_name": "Bombus species binomial epithet"})

    encoding = {
        name: {"zlib": True, "complevel": 4}
        for name in ds_out.data_vars
    }
    ds_out.to_netcdf(out_path, engine="netcdf4", encoding=encoding)
    print(f"  Saved {out_path}")

    # ============================================================
    # Option B: ALSO produce the nside=128 version (Tier-1 GLMM at
    # nside=64 unchanged; per-cell future predictors at full DestinE
    # nside=128 resolution). Per-species historical baselines are
    # inherited from the nside=64 parent (`parent_row_for_child`),
    # since the GLMM was calibrated per parent. Future TEI/PEI varies
    # per nside=128 child because climate is finer.
    # ============================================================
    print(f"\n  --- Option B: nside=128 climate predictors ---")
    meanT_future_128 = mmean128_c.mean(axis=0)               # (N_128,) degC
    annual_p_128 = np.empty((len(unique_years), N_128))
    for i, y in enumerate(unique_years):
        sel = years_p == y
        annual_p_128[i] = msum128[sel].sum(axis=0)
    annual_p_128_mm = annual_p_128 * 1000.0 * TP_HOURLY_TO_DAILY_FACTOR
    meanP_future_128 = annual_p_128_mm.mean(axis=0)          # (N_128,) mm/yr

    # Load Tier-1 nside=128 historical baseline (CRU TS sampled at
    # nside=128 cell centres) — the proper Option B reference.
    # Replaces the earlier parent-inheritance approach which inflated
    # η for cold-adapted species by mismatching the baseline scale.
    hist128_path = OUT_DIR / "climate_tei_pei_healpix_nside128.nc"
    if not hist128_path.exists():
        raise SystemExit(
            f"Missing {hist128_path}. Run healpix_port/04_climate_tei_pei_healpix.py "
            "to produce the Option-B nside=128 historical baseline."
        )
    hist128 = xr.open_dataset(hist128_path)
    TEI_bs_hist_128 = hist128["tei_bs"].values            # (n_spp, N_128) — proper baseline
    PEI_bs_hist_128 = hist128["pei_bs"].values
    avgtemp_bs_hist_128 = hist128["meanT_bs"].values      # (N_128,) — proper baseline
    avgprecip_bs_hist_128 = hist128["meanP_bs"].values

    with np.errstate(invalid="ignore", divide="ignore"):
        TEI_future_128 = (
            (meanT_future_128[np.newaxis, :] - T_min_spp[:, np.newaxis])
            / T_range[:, np.newaxis]
        )
        PEI_future_128 = (
            (meanP_future_128[np.newaxis, :] - P_min_spp[:, np.newaxis])
            / P_range[:, np.newaxis]
        )

    TEI_delta_128 = (TEI_future_128 - TEI_bs_hist_128).astype(np.float32)
    PEI_delta_128 = (PEI_future_128 - PEI_bs_hist_128).astype(np.float32)
    avgtemp_delta_128 = (meanT_future_128 - avgtemp_bs_hist_128).astype(np.float32)
    avgprecip_delta_128 = (meanP_future_128 - avgprecip_bs_hist_128).astype(np.float32)

    print(
        f"    meanT_future_128:  {meanT_future_128.min():.2f} .. "
        f"{meanT_future_128.max():.2f} degC  "
        f"(median {np.median(meanT_future_128):.2f}; "
        f"vs nside=64 median {np.median(meanT_future_64):.2f})"
    )
    print(
        f"    avgtemp_delta_128 median: "
        f"{np.median(avgtemp_delta_128):.3f} degC  "
        f"(vs nside=64 median {np.median(meanT_future_64 - avgtemp_bs_hist):.3f})"
    )

    out_path_128 = OUT_DIR / f"climate_tei_pei_future_{horizon}_healpix_nside128.nc"
    ds128 = xr.Dataset(
        data_vars={
            "tei_bs":          (("species", "cells"), TEI_bs_hist_128.astype(np.float32),
                                {"long_name": "baseline-period TEI (parent-inherited from Tier-1 nside=64)",
                                 "_FillValue": np.float32(np.nan)}),
            "tei_future":      (("species", "cells"), TEI_future_128.astype(np.float32),
                                {"long_name": "future-decade TEI (per-cell at nside=128)",
                                 "_FillValue": np.float32(np.nan)}),
            "tei_delta":       (("species", "cells"), TEI_delta_128,
                                {"long_name": "future − baseline TEI (per-cell at nside=128)",
                                 "_FillValue": np.float32(np.nan)}),
            "pei_bs":          (("species", "cells"), PEI_bs_hist_128.astype(np.float32),
                                {"long_name": "baseline-period PEI (parent-inherited from Tier-1 nside=64)",
                                 "_FillValue": np.float32(np.nan)}),
            "pei_future":      (("species", "cells"), PEI_future_128.astype(np.float32),
                                {"long_name": "future-decade PEI (per-cell at nside=128)",
                                 "_FillValue": np.float32(np.nan)}),
            "pei_delta":       (("species", "cells"), PEI_delta_128,
                                {"long_name": "future − baseline PEI (per-cell at nside=128)",
                                 "_FillValue": np.float32(np.nan)}),
            "meanT_future":    (("cells",), meanT_future_128.astype(np.float32),
                                {"long_name": "decade-mean monthly mean temperature",
                                 "units": "degC"}),
            "meanP_future":    (("cells",), meanP_future_128.astype(np.float32),
                                {"long_name": "decade-mean annual total precipitation",
                                 "units": "mm year-1"}),
            "avgtemp_bs":      (("cells",), avgtemp_bs_hist_128.astype(np.float32),
                                {"long_name": "baseline-period mean annual T (parent-inherited)",
                                 "units": "degC"}),
            "avgprecip_bs":    (("cells",), avgprecip_bs_hist_128.astype(np.float32),
                                {"long_name": "baseline-period annual total precipitation (parent-inherited)",
                                 "units": "mm year-1"}),
            "T_min_spp":       (("species",), T_min_spp.astype(np.float32),
                                {"long_name": "species cold thermal limit (Tier-1, unchanged)",
                                 "units": "degC"}),
            "T_max_spp":       (("species",), T_max_spp.astype(np.float32),
                                {"long_name": "species hot thermal limit (Tier-1, unchanged)",
                                 "units": "degC"}),
            "P_min_spp":       (("species",), P_min_spp.astype(np.float32),
                                {"long_name": "species dry precipitation limit (Tier-1, unchanged)",
                                 "units": "mm year-1"}),
            "P_max_spp":       (("species",), P_max_spp.astype(np.float32),
                                {"long_name": "species wet precipitation limit (Tier-1, unchanged)",
                                 "units": "mm year-1"}),
            "parent_row_in_nside64": (("cells",), parent_row_for_child.astype(np.int32),
                                {"long_name": "row index in nside=64 IBERIA_PIX_64 of this child's parent",
                                 "comment": "for sampling-effort & RE inheritance in 07"}),
        },
        coords={
            "species": np.array(species, dtype=object),
            "cell_ids": ("cells", IBERIA_PIX_128.astype(np.int64)),
        },
        attrs={
            "Conventions": "CF-1.10",
            "title": (
                f"DestinE Climate DT — Option B nside=128 future TEI/PEI "
                f"({decade_label}); baseline parent-inherited from Tier-1 nside=64"
            ),
            "history": (
                f"Created {date.today().isoformat()} by notebooks/06_destine_clean.py "
                "(Option B: GLMM at nside=64, climate predictors at nside=128)"
            ),
            "horizon": decade_label,
            "license_note": "DestinE Climate DT redistribution restricted; this file is gitignored.",
            **{**PROJECT_DGGS_ATTRS, "dggs_grid_refinement_level": 7},  # nside=128 → depth=7
            "tp_hourly_to_daily_factor": TP_HOURLY_TO_DAILY_FACTOR,
            "n_cells": int(N_128),
        },
    )
    ds128["cell_ids"].attrs.update({
        "long_name": "HEALPix NESTED pixel index (nside=128)",
    })
    ds128["species"].attrs.update({"long_name": "Bombus species binomial epithet"})
    enc128 = {name: {"zlib": True, "complevel": 4} for name in ds128.data_vars}
    ds128.to_netcdf(out_path_128, engine="netcdf4", encoding=enc128)
    print(f"    Saved {out_path_128}")

    # Diagnostics, restricted to cells with at least one species
    # historically observed (sampling_continent_healpix.nc active mask).
    sc = xr.open_dataset(OUT_DIR / "sampling_continent_healpix.nc")
    active = np.isfinite(sc["sampling_total"].values)
    print(
        f"  active cells (sampled at least once): {int(active.sum())}/"
        f"{N_64}"
    )
    print(
        f"  TEI_delta_future (species × active cells): "
        f"{np.nanmin(TEI_delta_future[:, active]):.3f} .. "
        f"{np.nanmax(TEI_delta_future[:, active]):.3f}  "
        f"(median {np.nanmedian(TEI_delta_future[:, active]):.3f})"
    )
    extreme = (TEI_future[:, active] > 1.0).sum()
    print(
        f"  cells × species exceeding historical hot-edge "
        f"(TEI_future > 1.0): {int(extreme):,}"
    )
