# HEALPix-NESTED nside=64 port of soroye_port/04_climate_tei_pei.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 04.

Computes per-species Thermal Exposure Index (TEI) and Precipitation
Exposure Index (PEI) on the Iberia HEALPix nside=64 NESTED substrate,
following the same Climatic Position Index (CPI) formula used by the
upstream / CEA pipeline:

    TEI[s, c] = (meanT[c] - T_min_spp[s]) / (T_max_spp[s] - T_min_spp[s])
    PEI[s, c] = (meanP[c] - P_min_spp[s]) / (P_max_spp[s] - P_min_spp[s])

where `meanT[c]` is the period-mean annual temperature at HEALPix cell
centre `c`, sampled from CRU TS 3.24.01 by **bilinear interpolation at
the cell centre** (single bilinear value -- no aggregation over CRU
pixels). The CRU TS grid (~0.5 deg) is much coarser than HEALPix
nside=64 (~92 km), so a centre-only sample is the natural and
defensible choice and matches what the upstream `04_climate_tei_pei.py`
does for CEA cell centres.

Per-species cold/hot/dry/wet limits are derived from the cells the
species occupied in the 1901-1974 baseline -- identical logic to
upstream; only the substrate (cell coordinates) changes.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / 'reference' / 'Bumblebee_repo_wbombusdat'
CLIM_DIR = REF / '0_ClimateData'
HEALPIX_PORT = ROOT / 'healpix_port'
OUT_DIR = HEALPIX_PORT / 'outputs_iberia'
PA_NPZ = OUT_DIR / 'presence_absence_healpix.npz'

# Period year sets (mirror upstream).
BASELINE_YEARS = list(range(1901, 1975))
RECENT_YEARS = list(range(2000, 2015))


# ---------------------------------------------------------------------------
# CRU TS loaders (verbatim from soroye_port/04_climate_tei_pei.py).

def load_cru_annual(var: str, agg: str) -> xr.DataArray:
    """Load + annually aggregate all CRU TS files for variable `var`.

    agg: 'mean' | 'sum' | 'min' | 'max'
    Returns DataArray (year, lat, lon) in degC (or mm for precip).
    """
    files = sorted(CLIM_DIR.glob(f'*{var}.dat.nc'))
    print(f'  {var}: {len(files)} decade files')
    pieces = []
    for f in files:
        ds = xr.open_dataset(f)
        v = ds[var]
        annual = getattr(v.resample(time='YE'), agg)()
        pieces.append(annual)
    annual_all = xr.concat(pieces, dim='time').sortby('time')
    annual_all = annual_all.assign_coords(time=annual_all['time.year']).rename({'time': 'year'})
    return annual_all


def bilinear_at_points(
    global_da: xr.DataArray,
    lat_pts: np.ndarray,
    lon_pts: np.ndarray,
) -> np.ndarray:
    """Bilinear-interpolate a global (lat, lon) DataArray at arbitrary
    (lat, lon) points. Returns a flat array aligned with the inputs.
    Identical algorithm to soroye_port/04_climate_tei_pei.py:bilinear_to_cea
    -- the only difference is the points come from HEALPix cell centres
    instead of CEA cell centres.
    """
    lat_coord = global_da['lat'].values
    lon_coord = global_da['lon'].values

    def to_idx(coord: np.ndarray, target: np.ndarray) -> np.ndarray:
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
    vals = map_coordinates(arr, [lat_idx, lon_idx], order=1,
                           mode='constant', cval=np.nan)
    return vals


# ---------------------------------------------------------------------------
# 1. Load Iberia cell-centre table from script 02 output.

print('Loading Iberia HEALPix cell-centre table ...')
pa = np.load(PA_NPZ, allow_pickle=True)
species = list(pa['species'])
prab_baseline = pa['prab_baseline']                 # (n_spp, n_cells)
iberia_lon = pa['iberia_lon'].astype(np.float64)    # cell-centre lon, [-180, 180]
iberia_lat = pa['iberia_lat'].astype(np.float64)    # cell-centre lat
n_cells = int(pa['n_cells'][0])
n_spp = len(species)
print(f'  {n_cells} cells, {n_spp} species')

# CRU TS uses lon in [-180, 180] (same as our wrapped HEALPix centres).

# ---------------------------------------------------------------------------
# 2. Load + annually aggregate CRU TS.

print('\nLoading CRU TS annual aggregates ...')
tmp_annual = load_cru_annual('tmp', 'mean')
pre_annual = load_cru_annual('pre', 'sum')
tmn_annual = load_cru_annual('tmn', 'min')
tmx_annual = load_cru_annual('tmx', 'max')
print(f'  Years available: {int(tmp_annual.year.min())}'
      f'..{int(tmp_annual.year.max())}')


# ---------------------------------------------------------------------------
# 3. Interpolate to HEALPix cell centres for each year.

def interp_years(da: xr.DataArray, years: list[int]) -> np.ndarray:
    """Return array shape (n_years, n_cells), bilinear at HEALPix
    cell centres."""
    out = np.full((len(years), n_cells), np.nan, dtype=np.float32)
    avail = set(int(y) for y in da.year.values)
    for i, yr in enumerate(years):
        if yr not in avail:
            continue
        layer = da.sel(year=yr)
        if 'latitude' in layer.dims:
            layer = layer.rename({'latitude': 'lat', 'longitude': 'lon'})
        out[i, :] = bilinear_at_points(layer, iberia_lat, iberia_lon)
    return out


print('\nInterpolating tmp (baseline + recent) to HEALPix cell centres ...')
tmp_bs_yr = interp_years(tmp_annual, BASELINE_YEARS)
tmp_rc_yr = interp_years(tmp_annual, RECENT_YEARS)
print(f'  baseline shape {tmp_bs_yr.shape}, recent shape {tmp_rc_yr.shape}')

print('Interpolating pre ...')
pre_bs_yr = interp_years(pre_annual, BASELINE_YEARS)
pre_rc_yr = interp_years(pre_annual, RECENT_YEARS)

print('Interpolating tmn / tmx (all years) for climate-cold/hot limits ...')
all_years = list(range(int(tmn_annual.year.min()),
                       int(tmn_annual.year.max()) + 1))
tmn_all_yr = interp_years(tmn_annual, all_years)
tmx_all_yr = interp_years(tmx_annual, all_years)

# ---------------------------------------------------------------------------
# 4. Period means (TEI uses CPI(mean_T) by linearity of CPI).

meanT_bs = np.nanmean(tmp_bs_yr, axis=0)
meanT_rc = np.nanmean(tmp_rc_yr, axis=0)
meanP_bs = np.nanmean(pre_bs_yr, axis=0)
meanP_rc = np.nanmean(pre_rc_yr, axis=0)

avgtemp_bs = meanT_bs
avgtemp_delta = meanT_rc - meanT_bs
avgprecip_bs = meanP_bs
avgprecip_delta = meanP_rc - meanP_bs

# ---------------------------------------------------------------------------
# 5. Per-cell climate cold/hot/dry/wet over all available years.

with np.errstate(invalid='ignore'):
    T_cold = np.nanmin(tmn_all_yr, axis=0)
    T_hot = np.nanmax(tmx_all_yr, axis=0)

pre_all_yr = interp_years(pre_annual, all_years)
with np.errstate(invalid='ignore'):
    P_dry = np.nanmin(pre_all_yr, axis=0)
    P_wet = np.nanmax(pre_all_yr, axis=0)

# ---------------------------------------------------------------------------
# 6. Per-species thermal + precip limits from baseline-occupied cells.

print('\nComputing per-species thermal + precip limits ...')
T_min_spp = np.full(n_spp, np.nan)
T_max_spp = np.full(n_spp, np.nan)
P_min_spp = np.full(n_spp, np.nan)
P_max_spp = np.full(n_spp, np.nan)
for s in range(n_spp):
    occupied = prab_baseline[s] == 1
    if not occupied.any():
        continue
    T_min_spp[s] = np.nanmin(T_cold[occupied])
    T_max_spp[s] = np.nanmax(T_hot[occupied])
    P_min_spp[s] = np.nanmin(P_dry[occupied])
    P_max_spp[s] = np.nanmax(P_wet[occupied])

# ---------------------------------------------------------------------------
# 7. TEI / PEI per species per cell per period.

print('Computing TEI / PEI ...')
T_range = T_max_spp - T_min_spp
P_range = P_max_spp - P_min_spp

TEI_bs = (meanT_bs[np.newaxis, :] - T_min_spp[:, np.newaxis]) / T_range[:, np.newaxis]
TEI_rc = (meanT_rc[np.newaxis, :] - T_min_spp[:, np.newaxis]) / T_range[:, np.newaxis]
TEI_delta = TEI_rc - TEI_bs

PEI_bs = (meanP_bs[np.newaxis, :] - P_min_spp[:, np.newaxis]) / P_range[:, np.newaxis]
PEI_rc = (meanP_rc[np.newaxis, :] - P_min_spp[:, np.newaxis]) / P_range[:, np.newaxis]
PEI_delta = PEI_rc - PEI_bs

# ---------------------------------------------------------------------------
# 8. Save.

np.savez_compressed(
    OUT_DIR / 'climate_tei_pei_healpix.npz',
    species=np.array(species),
    avgtemp_bs=avgtemp_bs.astype(np.float32),
    avgtemp_delta=avgtemp_delta.astype(np.float32),
    avgprecip_bs=avgprecip_bs.astype(np.float32),
    avgprecip_delta=avgprecip_delta.astype(np.float32),
    TEI_bs=TEI_bs.astype(np.float32),
    TEI_delta=TEI_delta.astype(np.float32),
    PEI_bs=PEI_bs.astype(np.float32),
    PEI_delta=PEI_delta.astype(np.float32),
    T_min_spp=T_min_spp.astype(np.float32),
    T_max_spp=T_max_spp.astype(np.float32),
    P_min_spp=P_min_spp.astype(np.float32),
    P_max_spp=P_max_spp.astype(np.float32),
    iberia_cells_hp=pa['iberia_cells_hp'].astype(np.uint64),
    n_cells=np.array([n_cells], dtype=np.int32),
)
print(f'\nSaved -> {OUT_DIR / "climate_tei_pei_healpix.npz"}')
print(f'\navgtemp_bs range: {np.nanmin(avgtemp_bs):.2f}..{np.nanmax(avgtemp_bs):.2f} degC')
print(f'avgtemp_delta range: {np.nanmin(avgtemp_delta):.2f}..{np.nanmax(avgtemp_delta):.2f} degC')
print(f'TEI_delta range: {np.nanmin(TEI_delta):.3f}..{np.nanmax(TEI_delta):.3f}')
print(f'Median TEI_delta at species-occupied baseline cells: '
      f'{np.nanmedian(TEI_delta[prab_baseline == 1]):.3f}')
