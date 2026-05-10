# HEALPix-NESTED nside=64 port of soroye_port/02_presence_absence.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 02.

Differs from the CEA port in one respect only: the underlying spatial
substrate. Instead of a 401x116 cylindrical-equal-area grid, we have a
flat list of Iberian HEALPix nside=64 NESTED cells (those whose centres
fall inside the lon -10..4, lat 35..44 bbox). Per-species presence,
inferred-absence, and species-richness logic is identical to the
upstream R / CEA port; only the indexing changes.

The presence/absence rule remains: in any (period, season) where any
species was observed at a cell, every other species without an
observation gets an explicit 0 (inferred absence).
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
HEALPIX_PORT = ROOT / 'healpix_port'
OUT_DIR = HEALPIX_PORT / 'outputs_iberia'
IN_CSV = OUT_DIR / 'bombus_clean_healpix.csv'

# HEALPix substrate constants (mirror script 01).
NSIDE = 64
DEPTH = 6
NPIX = 12 * NSIDE * NSIDE          # 49,152
ELLIPSOID = "WGS84"                # match DestinE Climate DT (per legacy-converters)

IBERIA_LON_MIN, IBERIA_LON_MAX = -10.0, 4.0
IBERIA_LAT_MIN, IBERIA_LAT_MAX = 35.0, 44.0


def _import_healpix_geo_nested():
    try:
        from healpix_geo import nested
        return nested
    except Exception:
        import healpix_geo as nested
        return nested


def pix_to_lonlat(ipix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pix -> (lon, lat) cell-center on WGS84 ellipsoid, lon wrapped to [-180, 180]."""
    nested = _import_healpix_geo_nested()
    ipix = np.asarray(ipix, dtype='uint64')
    lon, lat = nested.healpix_to_lonlat(ipix, DEPTH, ELLIPSOID)
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    return lon, lat


# ---------------------------------------------------------------------------
# 1. Build the flat Iberia cell list.
#
# We materialise the full Iberian cell set (NOT just cells with
# occurrences), so absences inferred from sampled-elsewhere logic stay
# meaningful at the substrate level. This mirrors the CEA pipeline,
# where the grid is the universe of cells, not the occurrence set.

all_pix = np.arange(NPIX, dtype='uint64')
all_lon, all_lat = pix_to_lonlat(all_pix)
iberia_mask = (
    (all_lon >= IBERIA_LON_MIN) & (all_lon <= IBERIA_LON_MAX)
    & (all_lat >= IBERIA_LAT_MIN) & (all_lat <= IBERIA_LAT_MAX)
)
iberia_cells_hp = all_pix[iberia_mask]                     # uint64 ipix
iberia_lon = all_lon[iberia_mask].astype(np.float64)
iberia_lat = all_lat[iberia_mask].astype(np.float64)
n_cells = len(iberia_cells_hp)
print(f'Iberia HEALPix nside={NSIDE} NESTED cells: {n_cells} '
      f'(bbox lon {IBERIA_LON_MIN}..{IBERIA_LON_MAX}, '
      f'lat {IBERIA_LAT_MIN}..{IBERIA_LAT_MAX})')

# Map: HEALPix ipix -> dense flat index in [0, n_cells)
ipix_to_idx = {int(p): i for i, p in enumerate(iberia_cells_hp)}


# ---------------------------------------------------------------------------
# 2. Load cleaned occurrences and assign each to its dense Iberia index.

print('\nLoading cleaned bombus data ...')
df = pd.read_csv(IN_CSV)
print(f'  {len(df):,} rows, {df["species"].nunique()} species')

df['cell_idx'] = df['cell_id_hp'].map(ipix_to_idx)
n_in = int(df['cell_idx'].notna().sum())
n_out = int(df['cell_idx'].isna().sum())
print(f'  {n_in:,} rows fall inside Iberia HEALPix mask, '
      f'{n_out:,} outside (rejected)')
df = df[df['cell_idx'].notna()].copy()
df['cell_idx'] = df['cell_idx'].astype(int)

species_list = sorted(df['species'].unique())
period_seasons = ['0_1', '0_2', '0_3', '3_1', '3_2', '3_3']


# ---------------------------------------------------------------------------
# 3. Per-season presence matrix [species x cell_idx], 1 / NaN.

print('\nBuilding presence matrices per (species x period_season) ...')

presence_key = (
    df.groupby(['period_season', 'species', 'cell_idx'])
      .size().reset_index(name='n_obs')
)

pre = {}
for ps in period_seasons:
    sub = presence_key[presence_key['period_season'] == ps]
    mat = np.full((len(species_list), n_cells), np.nan, dtype=np.float32)
    for _, row in sub.iterrows():
        spp_idx = species_list.index(row['species'])
        mat[spp_idx, int(row['cell_idx'])] = 1.0
    pre[ps] = mat
    print(f'  {ps}: {(mat == 1).sum():,} species x cell presences')


# ---------------------------------------------------------------------------
# 4. Species richness per period (for diagnostics, and for the npz).

print('\nComputing species richness per period ...')


def per_period_min(season_arrays: list[np.ndarray]) -> np.ndarray:
    stacked = np.stack(season_arrays, axis=0)
    with np.errstate(invalid='ignore'):
        any_pres = np.nanmin(stacked, axis=0)
    return any_pres


beedat_pr_baseline = per_period_min([pre['0_1'], pre['0_2'], pre['0_3']])
beedat_pr_recent = per_period_min([pre['3_1'], pre['3_2'], pre['3_3']])

sprich_baseline = np.nansum(beedat_pr_baseline, axis=0)
sprich_recent = np.nansum(beedat_pr_recent, axis=0)
sprich_baseline[sprich_baseline == 0] = np.nan
sprich_recent[sprich_recent == 0] = np.nan

print(f'  cells with any presence baseline: {np.isfinite(sprich_baseline).sum():,}')
print(f'  cells with any presence recent:   {np.isfinite(sprich_recent).sum():,}')


# ---------------------------------------------------------------------------
# 5. Inferred presence/absence: cells sampled anywhere (across all 6
#    period_seasons) get explicit 0 for species not observed there.

print('\nBuilding presence/absence (threshold = any species seen) ...')

total_sprich_cells = np.zeros(n_cells, dtype=float)
for ps in period_seasons:
    total_sprich_cells += np.nansum(pre[ps], axis=0)
sampled_anywhere = total_sprich_cells > 0

prab = {}
for ps in period_seasons:
    p = pre[ps].copy()
    mask = sampled_anywhere[np.newaxis, :] & np.isnan(p)
    p[mask] = 0.0
    prab[ps] = p
    print(f'  {ps}: {int((p == 0).sum()):,} inferred absences; '
          f'{int((p == 1).sum()):,} presences')


def per_period_max(prab_season_list):
    stacked = np.stack(prab_season_list, axis=0)
    with np.errstate(invalid='ignore'):
        m = np.nanmax(stacked, axis=0)
    return m


prab_baseline = per_period_max([prab['0_1'], prab['0_2'], prab['0_3']])
prab_recent = per_period_max([prab['3_1'], prab['3_2'], prab['3_3']])


# ---------------------------------------------------------------------------
# 6. Save as CF-compliant NetCDF.

# Stack 6 (period_season) x species x cell arrays from the per-season dicts.
presence_psc = np.stack(
    [pre[ps] for ps in period_seasons], axis=0
).astype(np.float32)                                       # (6, n_spp, n_cells)
prab_psc = np.stack(
    [prab[ps] for ps in period_seasons], axis=0
).astype(np.float32)

# Float fields where 0 means "absence inferred" must keep NaN as
# missing-data (CF-compliant via _FillValue) — float32 NaN already
# in arrays. Integer species_richness gets a sentinel for missing.
sprich_baseline_int = np.where(
    np.isfinite(sprich_baseline), sprich_baseline, -1,
).astype(np.int32)
sprich_recent_int = np.where(
    np.isfinite(sprich_recent), sprich_recent, -1,
).astype(np.int32)

ds = xr.Dataset(
    data_vars={
        'presence': (
            ('period_season', 'species', 'cell'),
            presence_psc,
            {
                'long_name': 'raw species presence per (period_season, species, cell)',
                'description': (
                    '1 = species observed in cell during period_season; '
                    'NaN = unknown (no observation). Inferred absences '
                    'are NOT applied here; see presence_absence.'
                ),
                '_FillValue': np.float32(np.nan),
            },
        ),
        'presence_absence': (
            ('period_season', 'species', 'cell'),
            prab_psc,
            {
                'long_name': 'inferred presence/absence per (period_season, species, cell)',
                'description': (
                    '1 = species observed in cell during period_season; '
                    '0 = inferred absence (some other Bombus seen in that '
                    'cell-period at any season); NaN = unknown'
                ),
                '_FillValue': np.float32(np.nan),
            },
        ),
        'prab_baseline': (
            ('species', 'cell'),
            prab_baseline.astype(np.float32),
            {
                'long_name': 'baseline-period (1901-1974) inferred presence/absence',
                'description': 'per-species per-cell max over baseline seasons',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'prab_recent': (
            ('species', 'cell'),
            prab_recent.astype(np.float32),
            {
                'long_name': 'recent-period (2000-2014) inferred presence/absence',
                'description': 'per-species per-cell max over recent seasons',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'species_richness_baseline': (
            ('cell',),
            sprich_baseline_int,
            {
                'long_name': 'number of species observed per cell in baseline period',
                'units': 'species',
                '_FillValue': np.int32(-1),
            },
        ),
        'species_richness_recent': (
            ('cell',),
            sprich_recent_int,
            {
                'long_name': 'number of species observed per cell in recent period',
                'units': 'species',
                '_FillValue': np.int32(-1),
            },
        ),
        'lon': (
            ('cell',),
            iberia_lon.astype(np.float32),
            {
                'long_name': 'HEALPix cell-centre longitude',
                'standard_name': 'longitude',
                'units': 'degrees_east',
            },
        ),
        'lat': (
            ('cell',),
            iberia_lat.astype(np.float32),
            {
                'long_name': 'HEALPix cell-centre latitude',
                'standard_name': 'latitude',
                'units': 'degrees_north',
            },
        ),
    },
    coords={
        'period_season': np.array(period_seasons, dtype='U5'),
        'species': np.array(species_list, dtype=object),
        'cell': iberia_cells_hp.astype(np.int64),
    },
    attrs={
        'Conventions': 'CF-1.10',
        'title': 'Iberian Bombus presence/absence on HEALPix nside=64 NESTED',
        'source': 'Soroye et al. 2020 method ported to HEALPix substrate',
        'history': (
            f'Created {date.today().isoformat()} by '
            'healpix_port/02_presence_absence_healpix.py'
        ),
        'healpix_nside': NSIDE,
        'healpix_depth': DEPTH,
        'healpix_scheme': 'NESTED',
        'healpix_lonlat_convention': 'WGS84, lon in [-180, 180]',
        'n_cells': n_cells,
    },
)

ds['cell'].attrs.update({
    'long_name': 'HEALPix NESTED pixel index (nside=64)',
    'description': 'Iberian-mask cell indices into the global HEALPix-NESTED nside=64 sphere',
})
ds['species'].attrs.update({'long_name': 'Bombus species binomial epithet'})
ds['period_season'].attrs.update({
    'long_name': 'period (0=baseline 1901-1974, 3=recent 2000-2014) and season (1-3)',
    'description': "Format '<period>_<season>'; periods are '0' (baseline) and '3' (recent)",
})

out_path = OUT_DIR / 'presence_absence_healpix.nc'
encoding = {
    'presence': {'zlib': True, 'complevel': 4},
    'presence_absence': {'zlib': True, 'complevel': 4},
    'prab_baseline': {'zlib': True, 'complevel': 4},
    'prab_recent': {'zlib': True, 'complevel': 4},
    'species_richness_baseline': {'zlib': True, 'complevel': 4},
    'species_richness_recent': {'zlib': True, 'complevel': 4},
    'lon': {'zlib': True, 'complevel': 4},
    'lat': {'zlib': True, 'complevel': 4},
}
ds.to_netcdf(out_path, engine='netcdf4', encoding=encoding)
print(f'\nSaved -> {out_path}')

# ---------------------------------------------------------------------------
# 7. Quick per-species summary.

print('\nPer-species cell counts (sample of 10 species):')
for i, spp in enumerate(species_list[:10]):
    n_bs = int((prab_baseline[i] == 1).sum())
    n_rc = int((prab_recent[i] == 1).sum())
    print(f'  {spp:<30} baseline={n_bs:4d}  recent={n_rc:4d}')
print(f'  ... ({len(species_list)} species total)')
