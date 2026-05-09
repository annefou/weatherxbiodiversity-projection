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
# # 06 — DestinE clean: future TEI / PEI on the Tier-1 CEA grid
#
# Aggregates the DestinE Climate DT daily extracts into the same CEA
# grid + same TEI / PEI definition as Tier 1 — so the GLMM coefficients
# fit in Tier 1 apply directly to the future predictors.
#
# **Critical:** species' thermal niche limits (`T_min_spp`, `T_max_spp`,
# `P_min_spp`, `P_max_spp`) are **historical** — derived in
# `soroye_port/04_climate_tei_pei.py` from the cells the species
# occupied in the 1901–1974 baseline. We hold them fixed; only the
# meanT / meanP per cell change with the future climate. This matches
# how `MapThermalStress.R` in the original Soroye codebase projects
# CPI under perturbed climate.
#
# Per Soroye 2020 / the upstream port, TEI is the **Climatic Position
# Index** evaluated at the period-mean annual temperature:
#
#     TEI[s, c] = (meanT[c] - T_min_spp[s]) / (T_max_spp[s] - T_min_spp[s])
#
# (Note: the user task brief described TEI as "fraction of months
# exceeding T_hot" — that is the alternative *exposure* definition.
# We follow the upstream port verbatim because that is what trained
# the GLMM coefficients we sample in 07.)
#
# Outputs (gitignored — see `.gitignore` Tier-2 section):
#
#   * `soroye_port/outputs_iberia/climate_tei_pei_future_2046_2055.npz`
#   * `soroye_port/outputs_iberia/climate_tei_pei_future_2076_2085.npz`

# %%
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import Transformer
from scipy.ndimage import map_coordinates

# %%
ROOT = Path("..").resolve()
DATA_DIR = ROOT / "data" / "destine"
PORT = ROOT / "soroye_port"
OUT_DIR = PORT / "outputs_iberia"

# %% [markdown]
# ## Re-use Tier-1 CEA grid construction
#
# Verbatim from `soroye_port/04_climate_tei_pei.py` so the output cell
# layout (401 × 116, 46 516 cells) matches Tier 1 byte-for-byte.

# %%
RES_M = 100_000
X_MIN, X_MAX = -20_037_507, 20_062_493
Y_MIN, Y_MAX = -5_263_885, 6_336_115

cea_to_ll = Transformer.from_crs(
    "+proj=cea +lat_ts=0 +lon_0=0 +ellps=WGS84",
    "EPSG:4326",
    always_xy=True,
)


def build_cea_cell_centers():
    n_x = (X_MAX - X_MIN) // RES_M
    n_y = (Y_MAX - Y_MIN) // RES_M
    x_centers = X_MIN + (np.arange(n_x) + 0.5) * RES_M
    y_centers = Y_MAX - (np.arange(n_y) + 0.5) * RES_M
    xx, yy = np.meshgrid(x_centers, y_centers)
    lon, lat = cea_to_ll.transform(xx.ravel(), yy.ravel())
    return n_x, n_y, lat.reshape((n_y, n_x)), lon.reshape((n_y, n_x))


n_x, n_y, cea_lat, cea_lon = build_cea_cell_centers()
n_cells = n_x * n_y
lat_flat = cea_lat.ravel()
lon_flat = cea_lon.ravel()
print(f"CEA grid: {n_x} × {n_y} = {n_cells:,} cells")


# %% [markdown]
# ## Load historical species niche limits (Tier 1) — held fixed

# %%
hist = np.load(OUT_DIR / "climate_tei_pei.npz", allow_pickle=True)
species = list(hist["species"])
n_spp = len(species)
T_min_spp = hist["T_min_spp"]
T_max_spp = hist["T_max_spp"]
P_min_spp = hist["P_min_spp"]
P_max_spp = hist["P_max_spp"]
TEI_bs_hist = hist["TEI_bs"]   # baseline-period TEI (n_spp, n_cells)
PEI_bs_hist = hist["PEI_bs"]
avgtemp_bs_hist = hist["avgtemp_bs"]
avgprecip_bs_hist = hist["avgprecip_bs"]
print(f"Historical limits loaded for {n_spp} species.")


# %% [markdown]
# ## Bilinear interpolation helper (lat/lon DataArray → CEA cells)
#
# Same shape and identical behaviour to the helper in
# `soroye_port/04_climate_tei_pei.py`.

# %%

def bilinear_to_cea(global_da: xr.DataArray, lat_pts: np.ndarray, lon_pts: np.ndarray) -> np.ndarray:
    """Bilinear-interpolate a global (lat, lon) DataArray onto the CEA cell centres."""
    if "latitude" in global_da.dims:
        global_da = global_da.rename({"latitude": "lat", "longitude": "lon"})
    lat_coord = global_da["lat"].values
    lon_coord = global_da["lon"].values

    def to_idx(coord, target):
        if coord[0] > coord[-1]:
            coord = coord[::-1]
            flipped = True
        else:
            flipped = False
        step = coord[1] - coord[0]
        idx = (target - coord[0]) / step
        if flipped:
            idx = (len(coord) - 1) - idx
        return idx

    lat_idx = to_idx(lat_coord, lat_pts)
    lon_idx = to_idx(lon_coord, lon_pts)
    arr = global_da.values.astype(float)
    return map_coordinates(arr, [lat_idx, lon_idx], order=1, mode="constant", cval=np.nan)


# %% [markdown]
# ## Per-horizon: load DestinE NetCDF → daily-to-monthly → period-mean → CEA → TEI/PEI

# %%
HORIZONS = ["2046_2055", "2076_2085"]


def _select_var(ds: xr.Dataset, candidates: list[str]) -> xr.DataArray:
    """Pick the first present variable from a candidate list. CHECK on first DestinE run."""
    for name in candidates:
        if name in ds.data_vars:
            return ds[name]
    raise KeyError(
        f"None of {candidates} in DestinE dataset. "
        f"CHECK: actual variable names are {list(ds.data_vars)}"
    )


def _ensure_celsius(da: xr.DataArray) -> xr.DataArray:
    """Convert K → °C if the variable looks like absolute temperature."""
    units = str(da.attrs.get("units", "")).lower()
    if units in ("k", "kelvin"):
        return da - 273.15
    if float(np.nanmean(da.values)) > 100.0:  # heuristic: temperatures in K range
        print(f"  CHECK: heuristic K→°C conversion applied for {da.name} (mean ~ {float(np.nanmean(da.values)):.1f})")
        return da - 273.15
    return da


def _ensure_mm_per_period(da: xr.DataArray) -> xr.DataArray:
    """Convert m → mm for total precipitation if needed."""
    units = str(da.attrs.get("units", "")).lower()
    if units in ("m", "metre", "meter"):
        return da * 1000.0
    return da


for horizon in HORIZONS:
    src = DATA_DIR / f"destine_iberia_{horizon}.nc"
    if not src.exists():
        print(f"[skip] {src.name} missing — run 05_destine_download.py first.")
        continue

    print(f"\n=== Horizon {horizon} — {src.name} ===")
    ds = xr.open_dataset(src)
    print(f"  variables: {list(ds.data_vars)}")
    print(f"  coords: {list(ds.coords)}")

    # CHECK: variable names depend on the DestinE catalogue. Common
    # GRIB-2 short names: '2t' / 'mx2t' / 'mn2t' / 'tp'. Try the most
    # plausible candidates and bail with a clear error otherwise.
    tmx = _ensure_celsius(_select_var(ds, ["mx2t", "tasmax", "tmax", "tx"]))
    tmn = _ensure_celsius(_select_var(ds, ["mn2t", "tasmin", "tmin", "tn"]))
    pre = _ensure_mm_per_period(_select_var(ds, ["tp", "pr", "precip", "pre"]))

    # Daily → monthly aggregates (max for tmax, min for tmin, sum for precip).
    # CHECK: time coord name; some DestinE outputs use "valid_time" or "step".
    time_dim = "time" if "time" in tmx.dims else list(tmx.dims)[0]
    print(f"  using time dim: {time_dim}")

    tmx_mon = tmx.resample({time_dim: "MS"}).max()
    tmn_mon = tmn.resample({time_dim: "MS"}).min()
    pre_mon = pre.resample({time_dim: "MS"}).sum()

    # Annual aggregates (matches Tier-1 contract: tmx→annual max,
    # tmn→annual min, pre→annual sum; tmp annual mean is built from
    # the monthly aggregates).
    tmx_ann = tmx_mon.resample({time_dim: "YS"}).max()
    tmn_ann = tmn_mon.resample({time_dim: "YS"}).min()
    pre_ann = pre_mon.resample({time_dim: "YS"}).sum()
    # Mean annual T from the monthly tmx/tmn midpoint (best surrogate
    # when DestinE returns only daily extremes; use mean(tmx,tmn) for
    # the period mean). CHECK: if DestinE also returns daily mean T
    # ('2t' / 'tas'), prefer that here.
    tmean_proxy_mon = (tmx_mon + tmn_mon) / 2.0
    tmean_ann = tmean_proxy_mon.resample({time_dim: "YS"}).mean()
    print("  CHECK: tmean computed as monthly midpoint of (tmx, tmn). "
          "If DestinE returns daily mean T directly, swap that in for "
          "exact comparability with CRU TS 'tmp'.")

    # Period mean over the decade — CPI uses the period-mean climate.
    meanT_future = tmean_ann.mean(dim=time_dim)
    meanP_future = pre_ann.mean(dim=time_dim)

    # Bilinear interpolation onto the CEA grid (same call signature
    # as Tier 1).
    print("  Interpolating onto CEA grid …")
    meanT_future_cea = bilinear_to_cea(meanT_future, lat_flat, lon_flat).astype(np.float32)
    meanP_future_cea = bilinear_to_cea(meanP_future, lat_flat, lon_flat).astype(np.float32)

    # Future-decade TEI / PEI per (species × cell), reusing historical
    # niche limits.
    T_range = T_max_spp - T_min_spp
    P_range = P_max_spp - P_min_spp
    with np.errstate(invalid="ignore", divide="ignore"):
        TEI_future = (
            (meanT_future_cea[np.newaxis, :] - T_min_spp[:, np.newaxis])
            / T_range[:, np.newaxis]
        )
        PEI_future = (
            (meanP_future_cea[np.newaxis, :] - P_min_spp[:, np.newaxis])
            / P_range[:, np.newaxis]
        )

    # Deltas relative to the historical baseline (so they line up with
    # the GLMM's `sc_TEI_delta` / `sc_PEI_delta` predictor — both are
    # already-historical-anchored).
    TEI_delta_future = (TEI_future - TEI_bs_hist).astype(np.float32)
    PEI_delta_future = (PEI_future - PEI_bs_hist).astype(np.float32)

    avgtemp_delta_future = (meanT_future_cea - avgtemp_bs_hist).astype(np.float32)
    avgprecip_delta_future = (meanP_future_cea - avgprecip_bs_hist).astype(np.float32)

    # Save with same key shape as Tier-1 npz so consumers can swap dictionaries.
    out_path = OUT_DIR / f"climate_tei_pei_future_{horizon}.npz"
    np.savez_compressed(
        out_path,
        species=np.array(species),
        # Baseline TEI/PEI carried over verbatim — projection step uses
        # the same `sc_TEI_bs` / `sc_PEI_bs` as Tier 1, so the predictor
        # row aligns at the species' historical-context baseline.
        TEI_bs=TEI_bs_hist.astype(np.float32),
        PEI_bs=PEI_bs_hist.astype(np.float32),
        TEI_delta=TEI_delta_future,
        PEI_delta=PEI_delta_future,
        TEI_future=TEI_future.astype(np.float32),
        PEI_future=PEI_future.astype(np.float32),
        avgtemp_bs=avgtemp_bs_hist,
        avgtemp_delta=avgtemp_delta_future,
        avgprecip_bs=avgprecip_bs_hist,
        avgprecip_delta=avgprecip_delta_future,
        meanT_future=meanT_future_cea,
        meanP_future=meanP_future_cea,
        T_min_spp=T_min_spp,
        T_max_spp=T_max_spp,
        P_min_spp=P_min_spp,
        P_max_spp=P_max_spp,
        n_x=n_x,
        n_y=n_y,
        horizon=horizon,
    )
    print(f"  Saved {out_path}")

    # Diagnostics — aggregate over the active Iberian cells only.
    sc = np.load(OUT_DIR / "sampling_continent.npz", allow_pickle=True)
    active = ~np.isnan(sc["samp_total"])
    print(f"  active cells: {active.sum()}")
    print(
        f"  meanT_future (active cells): "
        f"{np.nanmin(meanT_future_cea[active]):.2f} .. "
        f"{np.nanmax(meanT_future_cea[active]):.2f} °C  "
        f"(median {np.nanmedian(meanT_future_cea[active]):.2f})"
    )
    print(
        f"  avgtemp_delta vs baseline (active cells): "
        f"{np.nanmin(avgtemp_delta_future[active]):.2f} .. "
        f"{np.nanmax(avgtemp_delta_future[active]):.2f} °C  "
        f"(median {np.nanmedian(avgtemp_delta_future[active]):.2f})"
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
