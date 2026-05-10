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

from datetime import date
from pathlib import Path

import numpy as np
import xarray as xr

from _dggs_metadata import PROJECT_DGGS_ATTRS
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / 'reference' / 'Bumblebee_repo_wbombusdat'
CLIM_DIR = REF / '0_ClimateData'
HEALPIX_PORT = ROOT / 'healpix_port'
OUT_DIR = HEALPIX_PORT / 'outputs_iberia'
PA_NC = OUT_DIR / 'presence_absence_healpix.nc'

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
pa = xr.open_dataset(PA_NC)
species = [str(s) for s in pa['species'].values]
prab_baseline = pa['prab_baseline'].values         # (n_spp, n_cells)
iberia_lon = pa['lon'].values.astype(np.float64)   # cell-centre lon, [-180, 180]
iberia_lat = pa['lat'].values.astype(np.float64)   # cell-centre lat
iberia_cells_hp = pa['cell_ids'].values.astype(np.uint64)
n_cells = int(pa.sizes['cells'])
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

TEI_rc_arr = (TEI_bs + TEI_delta).astype(np.float32)
PEI_rc_arr = (PEI_bs + PEI_delta).astype(np.float32)
meanT_rc_arr = meanT_rc.astype(np.float32)
meanP_rc_arr = meanP_rc.astype(np.float32)

ds = xr.Dataset(
    data_vars={
        'tei_bs': (
            ('species', 'cells'),
            TEI_bs.astype(np.float32),
            {
                'long_name': 'Climatic Position Index (thermal) baseline 1901-1974',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'tei_rc': (
            ('species', 'cells'),
            TEI_rc_arr,
            {
                'long_name': 'Climatic Position Index (thermal) recent 2000-2014',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'tei_delta': (
            ('species', 'cells'),
            TEI_delta.astype(np.float32),
            {
                'long_name': 'Delta thermal CPI = recent minus baseline',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'pei_bs': (
            ('species', 'cells'),
            PEI_bs.astype(np.float32),
            {
                'long_name': 'Climatic Position Index (precipitation) baseline 1901-1974',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'pei_rc': (
            ('species', 'cells'),
            PEI_rc_arr,
            {
                'long_name': 'Climatic Position Index (precipitation) recent 2000-2014',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'pei_delta': (
            ('species', 'cells'),
            PEI_delta.astype(np.float32),
            {
                'long_name': 'Delta precipitation CPI = recent minus baseline',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'T_min_spp': (
            ('species',),
            T_min_spp.astype(np.float32),
            {
                'long_name': 'species-specific cold thermal limit (5-lowest baseline monthly Tmin avg)',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'T_max_spp': (
            ('species',),
            T_max_spp.astype(np.float32),
            {
                'long_name': 'species-specific hot thermal limit (max of baseline monthly Tmax)',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'P_min_spp': (
            ('species',),
            P_min_spp.astype(np.float32),
            {
                'long_name': 'species-specific dry-period precipitation limit',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'P_max_spp': (
            ('species',),
            P_max_spp.astype(np.float32),
            {
                'long_name': 'species-specific wet-period precipitation limit',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'meanT_bs': (
            ('cells',),
            meanT_bs.astype(np.float32),
            {
                'long_name': 'baseline-period mean annual temperature',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'meanT_rc': (
            ('cells',),
            meanT_rc_arr,
            {
                'long_name': 'recent-period mean annual temperature',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'meanP_bs': (
            ('cells',),
            meanP_bs.astype(np.float32),
            {
                'long_name': 'baseline-period mean annual total precipitation',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'meanP_rc': (
            ('cells',),
            meanP_rc_arr,
            {
                'long_name': 'recent-period mean annual total precipitation',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
        # Aliased convenience names (matching the legacy npz keys for
        # downstream scripts that still refer to avgtemp_*).
        'avgtemp_bs': (
            ('cells',),
            avgtemp_bs.astype(np.float32),
            {
                'long_name': 'baseline mean annual temperature (alias of meanT_bs)',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'avgtemp_delta': (
            ('cells',),
            avgtemp_delta.astype(np.float32),
            {
                'long_name': 'delta annual mean temperature (recent minus baseline)',
                'units': 'degC',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'avgprecip_bs': (
            ('cells',),
            avgprecip_bs.astype(np.float32),
            {
                'long_name': 'baseline mean annual total precipitation (alias of meanP_bs)',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'avgprecip_delta': (
            ('cells',),
            avgprecip_delta.astype(np.float32),
            {
                'long_name': 'delta annual total precipitation (recent minus baseline)',
                'units': 'mm',
                '_FillValue': np.float32(np.nan),
            },
        ),
    },
    coords={
        'species': np.array(species, dtype=object),
        "cell_ids": ("cells", iberia_cells_hp.astype(np.int64)),    },
    attrs={
        'Conventions': 'CF-1.10',
        'title': (
            'Climatic Position Index (TEI / PEI) per species per cell on '
            'Iberian HEALPix nside=64 NESTED'
        ),
        'source': 'Soroye et al. 2020 method ported to HEALPix substrate',
        'crs_source': 'CRU TS 3.24.01 (Soroye Figshare)',
        'history': (
            f'Created {date.today().isoformat()} by '
            'healpix_port/04_climate_tei_pei_healpix.py'
        ),
        **PROJECT_DGGS_ATTRS,    # DGGS Zarr Convention v1 — see _dggs_metadata.py
        'n_cells': n_cells,
        'baseline_period': '1901-1974',
        'recent_period': '2000-2014',
    },
)
ds['cell_ids'].attrs.update({
    'long_name': 'HEALPix NESTED pixel index (nside=64)',
})
ds['species'].attrs.update({'long_name': 'Bombus species binomial epithet'})

out_path = OUT_DIR / 'climate_tei_pei_healpix.nc'
encoding = {
    name: {'zlib': True, 'complevel': 4}
    for name in ds.data_vars
}
ds.to_netcdf(out_path, engine='netcdf4', encoding=encoding)
print(f'\nSaved -> {out_path}')


# =====================================================================
# OPTION B (proper) — nside=128 historical baseline
# =====================================================================
# Sample CRU TS at nside=128 cell centres (440 Iberian cells), compute
# the same Climatic Position Index machinery, derive standardisation
# constants from the nside=128 distribution. The species' thermal +
# precipitation niche limits (T_min_spp, T_max_spp, P_min_spp, P_max_spp)
# are SPECIES-LEVEL — kept identical to the nside=64 fit so the GLMM's
# coefficient β scaling is preserved at species granularity. What
# CHANGES at nside=128 is the per-cell baseline climate (sampled from
# CRU at nside=128 cell centres, NOT parent-inherited from nside=64).
#
# This file is consumed by notebooks/06_destine_clean.py to compute
# proper TEI_delta = TEI_future_128 - TEI_bs_128 (both at nside=128),
# and by notebooks/07_projection.py to z-score future predictors with
# the nside=128 standardisation reference.

print('\n--- Option B extension: nside=128 baseline climate ---')

PRECOMP = ROOT / 'data' / 'precomputed'
iberia_pix_128 = np.load(PRECOMP / 'iberia_pix_nside128_nested.npy').astype(np.uint64)
n_cells_128 = len(iberia_pix_128)

import healpix_geo  # noqa: E402
lon_128_arr, lat_128_arr = healpix_geo.nested.healpix_to_lonlat(
    iberia_pix_128, 7, 'WGS84',
)
lon_128_arr = np.where(lon_128_arr > 180.0, lon_128_arr - 360.0, lon_128_arr)
print(f'  {n_cells_128} nside=128 cell centres computed (WGS84)')


def interp_years_128(da: xr.DataArray, years: list[int]) -> np.ndarray:
    out = np.full((len(years), n_cells_128), np.nan, dtype=np.float32)
    avail = set(int(y) for y in da.year.values)
    for i, yr in enumerate(years):
        if yr not in avail:
            continue
        layer = da.sel(year=yr)
        if 'latitude' in layer.dims:
            layer = layer.rename({'latitude': 'lat', 'longitude': 'lon'})
        out[i, :] = bilinear_at_points(layer, lat_128_arr, lon_128_arr)
    return out


print('  Interpolating CRU TS to nside=128 cell centres ...')
tmp_bs_128 = interp_years_128(tmp_annual, BASELINE_YEARS)
tmp_rc_128 = interp_years_128(tmp_annual, RECENT_YEARS)
pre_bs_128 = interp_years_128(pre_annual, BASELINE_YEARS)
pre_rc_128 = interp_years_128(pre_annual, RECENT_YEARS)

meanT_bs_128 = np.nanmean(tmp_bs_128, axis=0)
meanT_rc_128 = np.nanmean(tmp_rc_128, axis=0)
meanP_bs_128 = np.nanmean(pre_bs_128, axis=0)
meanP_rc_128 = np.nanmean(pre_rc_128, axis=0)

# Per-species TEI/PEI at nside=128 (species niche limits unchanged from nside=64).
TEI_bs_128 = (meanT_bs_128[np.newaxis, :] - T_min_spp[:, np.newaxis]) / T_range[:, np.newaxis]
TEI_rc_128 = (meanT_rc_128[np.newaxis, :] - T_min_spp[:, np.newaxis]) / T_range[:, np.newaxis]
TEI_delta_128 = TEI_rc_128 - TEI_bs_128

PEI_bs_128 = (meanP_bs_128[np.newaxis, :] - P_min_spp[:, np.newaxis]) / P_range[:, np.newaxis]
PEI_rc_128 = (meanP_rc_128[np.newaxis, :] - P_min_spp[:, np.newaxis]) / P_range[:, np.newaxis]
PEI_delta_128 = PEI_rc_128 - PEI_bs_128

# Build per-species observation mask at nside=128 via parent-inheritance
# (each child gets its parent's historical observation status — there
# is no nside=128 GBIF cleaning, so this is the best available signal
# for "which species was historically present in this cell").
parent_to_row64 = {int(p): i for i, p in enumerate(iberia_cells_hp.astype(np.int64))}
parents_of_children = (iberia_pix_128.astype(np.int64) >> 2)
parent_row_128 = np.array(
    [parent_to_row64[int(p)] for p in parents_of_children], dtype=np.int64,
)
prab_baseline_128 = prab_baseline[:, parent_row_128]    # (n_spp, 440)


def std_over_observed(arr_128: np.ndarray, prab_128: np.ndarray) -> tuple[float, float]:
    """μ + ddof=1 σ over (species × cell) entries where the species was
    historically observed in the parent. Matches the Tier-1 standardisation
    reference shape (species × cell rows, observed-only)."""
    mask = (prab_128 > 0) & np.isfinite(arr_128)
    vals = arr_128[mask]
    if vals.size < 2:
        return float('nan'), float('nan')
    return float(vals.mean()), float(vals.std(ddof=1))


mu_TEI_bs_128, sd_TEI_bs_128 = std_over_observed(TEI_bs_128, prab_baseline_128)
mu_TEI_delta_128, sd_TEI_delta_128 = std_over_observed(TEI_delta_128, prab_baseline_128)
mu_PEI_bs_128, sd_PEI_bs_128 = std_over_observed(PEI_bs_128, prab_baseline_128)
mu_PEI_delta_128, sd_PEI_delta_128 = std_over_observed(PEI_delta_128, prab_baseline_128)

print('  Standardisation constants (nside=128 vs nside=64 in attrs):')
print(f'    TEI_bs:    mu={mu_TEI_bs_128:.4f}    sd={sd_TEI_bs_128:.4f}')
print(f'    TEI_delta: mu={mu_TEI_delta_128:.4f} sd={sd_TEI_delta_128:.4f}')
print(f'    PEI_bs:    mu={mu_PEI_bs_128:.4f}    sd={sd_PEI_bs_128:.4f}')
print(f'    PEI_delta: mu={mu_PEI_delta_128:.4f} sd={sd_PEI_delta_128:.4f}')

ds128 = xr.Dataset(
    data_vars={
        'tei_bs': (('species', 'cells'), TEI_bs_128.astype(np.float32),
                   {'long_name': 'CPI thermal baseline at nside=128 (CRU-sampled at cell centres)',
                    '_FillValue': np.float32(np.nan)}),
        'tei_rc': (('species', 'cells'), TEI_rc_128.astype(np.float32),
                   {'long_name': 'CPI thermal recent at nside=128',
                    '_FillValue': np.float32(np.nan)}),
        'tei_delta': (('species', 'cells'), TEI_delta_128.astype(np.float32),
                      {'long_name': 'CPI thermal delta = recent − baseline at nside=128',
                       '_FillValue': np.float32(np.nan)}),
        'pei_bs': (('species', 'cells'), PEI_bs_128.astype(np.float32),
                   {'long_name': 'CPI precipitation baseline at nside=128',
                    '_FillValue': np.float32(np.nan)}),
        'pei_rc': (('species', 'cells'), PEI_rc_128.astype(np.float32),
                   {'long_name': 'CPI precipitation recent at nside=128',
                    '_FillValue': np.float32(np.nan)}),
        'pei_delta': (('species', 'cells'), PEI_delta_128.astype(np.float32),
                      {'long_name': 'CPI precipitation delta at nside=128',
                       '_FillValue': np.float32(np.nan)}),
        'meanT_bs': (('cells',), meanT_bs_128.astype(np.float32),
                     {'long_name': 'baseline mean annual temperature at nside=128 cell centre',
                      'units': 'degC'}),
        'meanT_rc': (('cells',), meanT_rc_128.astype(np.float32),
                     {'long_name': 'recent mean annual temperature at nside=128 cell centre',
                      'units': 'degC'}),
        'meanP_bs': (('cells',), meanP_bs_128.astype(np.float32),
                     {'long_name': 'baseline annual total precipitation at nside=128 cell centre',
                      'units': 'mm year-1'}),
        'meanP_rc': (('cells',), meanP_rc_128.astype(np.float32),
                     {'long_name': 'recent annual total precipitation at nside=128 cell centre',
                      'units': 'mm year-1'}),
        'lon': (('cells',), lon_128_arr.astype(np.float32),
                {'standard_name': 'longitude', 'units': 'degrees_east'}),
        'lat': (('cells',), lat_128_arr.astype(np.float32),
                {'standard_name': 'latitude', 'units': 'degrees_north'}),
        'parent_row_in_nside64': (('cells',), parent_row_128.astype(np.int32),
                                  {'long_name': 'row index in nside=64 IBERIA_PIX_64 of this child\'s parent'}),
        'T_min_spp': (('species',), T_min_spp.astype(np.float32),
                      {'long_name': 'species cold thermal limit (Tier-1 nside=64-derived; species-level, unchanged)'}),
        'T_max_spp': (('species',), T_max_spp.astype(np.float32),
                      {'long_name': 'species hot thermal limit (species-level, unchanged)'}),
        'P_min_spp': (('species',), P_min_spp.astype(np.float32),
                      {'long_name': 'species dry precipitation limit (species-level, unchanged)'}),
        'P_max_spp': (('species',), P_max_spp.astype(np.float32),
                      {'long_name': 'species wet precipitation limit (species-level, unchanged)'}),
    },
    coords={
        'species': np.array(species, dtype=object),
        'cell_ids': ('cells', iberia_pix_128.astype(np.int64)),
    },
    attrs={
        'Conventions': 'CF-1.10',
        'title': 'Tier-1 historical CPI at HEALPix nside=128 NESTED — Option B baseline',
        'source': 'CRU TS 3.24.01 sampled at nside=128 cell centres (WGS84)',
        'history': f'Created {date.today().isoformat()} by healpix_port/04_climate_tei_pei_healpix.py',
        **{**PROJECT_DGGS_ATTRS, 'dggs_grid_refinement_level': 7},
        'n_cells': int(n_cells_128),
        'baseline_period': '1901-1974',
        'recent_period': '2000-2014',
        # Standardisation constants for Option B z-scoring at nside=128.
        # Used by notebooks/07_projection.py when applying Tier-1 GLMM
        # coefficients to nside=128 predictors.
        'std_TEI_bs_mu':    mu_TEI_bs_128,
        'std_TEI_bs_sd':    sd_TEI_bs_128,
        'std_TEI_delta_mu': mu_TEI_delta_128,
        'std_TEI_delta_sd': sd_TEI_delta_128,
        'std_PEI_bs_mu':    mu_PEI_bs_128,
        'std_PEI_bs_sd':    sd_PEI_bs_128,
        'std_PEI_delta_mu': mu_PEI_delta_128,
        'std_PEI_delta_sd': sd_PEI_delta_128,
        'std_reference': (
            'Computed at nside=128 over (species × cell) entries where the '
            'species was historically observed at the parent nside=64 cell. '
            'Different from Tier-1 nside=64 reference; applying Tier-1 GLMM '
            'coefficients to these z-scores assumes the GLMM\'s coefficient '
            'magnitude is comparable across z-score scales when applied at '
            'matched substrate scale.'
        ),
    },
)
ds128['cell_ids'].attrs.update({
    'long_name': 'HEALPix NESTED pixel index (nside=128)',
})
ds128['species'].attrs.update({'long_name': 'Bombus species binomial epithet'})

out_path_128 = OUT_DIR / 'climate_tei_pei_healpix_nside128.nc'
enc128 = {name: {'zlib': True, 'complevel': 4} for name in ds128.data_vars}
ds128.to_netcdf(out_path_128, engine='netcdf4', encoding=enc128)
print(f'  Saved -> {out_path_128}')
print(f'\navgtemp_bs range: {np.nanmin(avgtemp_bs):.2f}..{np.nanmax(avgtemp_bs):.2f} degC')
print(f'avgtemp_delta range: {np.nanmin(avgtemp_delta):.2f}..{np.nanmax(avgtemp_delta):.2f} degC')
print(f'TEI_delta range: {np.nanmin(TEI_delta):.3f}..{np.nanmax(TEI_delta):.3f}')
print(f'Median TEI_delta at species-occupied baseline cells: '
      f'{np.nanmedian(TEI_delta[prab_baseline == 1]):.3f}')
