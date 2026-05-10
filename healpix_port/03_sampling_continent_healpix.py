# HEALPix-NESTED nside=64 port of soroye_port/03_sampling_continent.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 03.

Computes per-cell sampling effort (count of distinct LYIDs per cell per
period_season, summed across the six season rasters) and the continent
code. Continent is constant = 2 (Iberia/Europe), so the continent
"raster" is trivial here, but we keep the field so the downstream
regression script can use the same code path as the upstream / CEA
version.
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
PA_NC = OUT_DIR / 'presence_absence_healpix.nc'

# ---------------------------------------------------------------------------
# 1. Load cleaned data and the Iberia cell list (the universe of sites).

print('Loading cleaned data + Iberia cell list ...')
df = pd.read_csv(IN_CSV)

pa = xr.open_dataset(PA_NC)
iberia_cells_hp = pa['cell'].values.astype(np.uint64)
n_cells = int(pa.sizes['cell'])
print(f'  Iberia cells: {n_cells}, occurrences: {len(df):,}')

# Map ipix -> dense Iberia index, mirroring script 02.
ipix_to_idx = {int(p): i for i, p in enumerate(iberia_cells_hp)}
df['cell_idx'] = df['cell_id_hp'].map(ipix_to_idx)
df = df[df['cell_idx'].notna()].copy()
df['cell_idx'] = df['cell_idx'].astype(int)

# ---------------------------------------------------------------------------
# 2. Unique LYID rows per (continent, lon, lat, period_season, species, LYID)
#    -- mirror upstream R script 3 line 27.

lyid_df = (
    df.groupby(
        ['continent', 'longitude', 'latitude', 'period_season',
         'species', 'LYID', 'cell_idx'],
        as_index=False,
    ).size().rename(columns={'size': 'nobs'})
)
print(f'  unique LYID rows: {len(lyid_df):,}')

# ---------------------------------------------------------------------------
# 3. Per (period_season, cell_idx) count of distinct LYIDs.

print('\nComputing sampling per season ...')
period_seasons = ['0_1', '0_2', '0_3', '3_1', '3_2', '3_3']

samp_seasons: dict[str, np.ndarray] = {}
for ps in period_seasons:
    sub = lyid_df[lyid_df['period_season'] == ps]
    count = np.full(n_cells, np.nan, dtype=np.float32)
    if len(sub) > 0:
        cell_counts = sub.groupby('cell_idx').size()
        count[cell_counts.index.values] = cell_counts.values.astype(np.float32)
    samp_seasons[ps] = count

# Cells sampled in any season -> fill NaN-in-this-season with 0
all_cells_counts = np.full(n_cells, np.nan, dtype=np.float32)
grouped_all = lyid_df.groupby('cell_idx').size()
all_cells_counts[grouped_all.index.values] = grouped_all.values.astype(np.float32)

for ps in period_seasons:
    v = samp_seasons[ps]
    fill_mask = (all_cells_counts > 0) & np.isnan(v)
    v[fill_mask] = 0.0
    samp_seasons[ps] = v
    print(f'  {ps}: mean={np.nanmean(v):.2f}  '
          f'nonzero={int((v > 0).sum()):,}  zeros={int((v == 0).sum()):,}')

# ---------------------------------------------------------------------------
# 4. Sum across 3 seasons in each period.

samp_baseline = np.nansum(
    np.stack([samp_seasons['0_1'], samp_seasons['0_2'], samp_seasons['0_3']], axis=0),
    axis=0,
)
samp_recent = np.nansum(
    np.stack([samp_seasons['3_1'], samp_seasons['3_2'], samp_seasons['3_3']], axis=0),
    axis=0,
)

any_bs = ~np.isnan(np.stack([samp_seasons['0_1'], samp_seasons['0_2'], samp_seasons['0_3']])).all(axis=0)
any_rc = ~np.isnan(np.stack([samp_seasons['3_1'], samp_seasons['3_2'], samp_seasons['3_3']])).all(axis=0)
samp_baseline[~any_bs] = np.nan
samp_recent[~any_rc] = np.nan

print(f'\nBaseline sampling: {int(np.isfinite(samp_baseline).sum()):,} cells, '
      f'total LYIDs {int(np.nansum(samp_baseline)):,}')
print(f'Recent sampling:   {int(np.isfinite(samp_recent).sum()):,} cells, '
      f'total LYIDs {int(np.nansum(samp_recent)):,}')

# ---------------------------------------------------------------------------
# 5. Continent raster -- constant 2 (Europe) for Iberian cells with sampling.

print('\nBuilding continent raster ...')
continent = np.full(n_cells, np.nan, dtype=np.float32)
agg = lyid_df.groupby('cell_idx')['continent'].mean()
continent[agg.index.values] = agg.values.astype(np.float32)
print(f'  cells with continent assigned: {int(np.isfinite(continent).sum()):,}')

# ---------------------------------------------------------------------------
# 6. Total sampling = sum of all 6 seasons (matches upstream R `sampling`).

samp_stack = np.stack([samp_seasons[ps] for ps in period_seasons], axis=0)
samp_total = np.nansum(samp_stack, axis=0)
samp_total[samp_total == 0] = np.nan
print(f'Total sampling (all 6 seasons): cells with sampling = '
      f'{int(np.isfinite(samp_total).sum()):,}, '
      f'total LYIDs = {int(np.nansum(samp_total)):,}')

# Continent: int8 with -1 sentinel for missing (NaN in float continent).
continent_int8 = np.where(
    np.isfinite(continent), continent, -1,
).astype(np.int8)

ds = xr.Dataset(
    data_vars={
        'sampling': (
            ('period_season', 'cell'),
            samp_stack.astype(np.float32),
            {
                'long_name': 'number of distinct LYID-equivalents per cell per period_season',
                'units': '1',
                'description': (
                    'Per (period, season) count of unique '
                    '(continent, lon, lat, species, LYID) combinations '
                    'in the cleaned Bombus occurrence table. '
                    'NaN where the cell was never sampled in any of '
                    'the six period_seasons.'
                ),
                '_FillValue': np.float32(np.nan),
            },
        ),
        'sampling_baseline': (
            ('cell',),
            samp_baseline.astype(np.float32),
            {
                'long_name': 'sampling effort summed across baseline-period seasons',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'sampling_recent': (
            ('cell',),
            samp_recent.astype(np.float32),
            {
                'long_name': 'sampling effort summed across recent-period seasons',
                'units': '1',
                '_FillValue': np.float32(np.nan),
            },
        ),
        'sampling_total': (
            ('cell',),
            samp_total.astype(np.float32),
            {
                'long_name': 'total sampling effort summed over all six period_seasons',
                'units': '1',
                'description': (
                    'Per-cell total LYID-count across all six (period, season) '
                    'combinations; matches upstream R `sampling`.'
                ),
                '_FillValue': np.float32(np.nan),
            },
        ),
        'continent': (
            ('cell',),
            continent_int8,
            {
                'long_name': 'continent index per cell (1=North America, 2=Europe)',
                'flag_values': np.array([1, 2], dtype=np.int8),
                'flag_meanings': 'north_america europe',
                '_FillValue': np.int8(-1),
            },
        ),
    },
    coords={
        'period_season': np.array(period_seasons, dtype='U5'),
        'cell': iberia_cells_hp.astype(np.int64),
    },
    attrs={
        'Conventions': 'CF-1.10',
        'title': 'Per-cell sampling effort + continent on Iberian HEALPix nside=64 NESTED',
        'source': 'Soroye et al. 2020 method ported to HEALPix substrate',
        'history': (
            f'Created {date.today().isoformat()} by '
            'healpix_port/03_sampling_continent_healpix.py'
        ),
        'healpix_nside': 64,
        'healpix_depth': 6,
        'healpix_scheme': 'NESTED',
        'healpix_lonlat_convention': 'WGS84, lon in [-180, 180]',
        'n_cells': n_cells,
    },
)
ds['cell'].attrs.update({
    'long_name': 'HEALPix NESTED pixel index (nside=64)',
})
ds['period_season'].attrs.update({
    'long_name': 'period_season label',
    'description': "Format '<period>_<season>'; periods 0=baseline, 3=recent.",
})

out_path = OUT_DIR / 'sampling_continent_healpix.nc'
encoding = {
    'sampling': {'zlib': True, 'complevel': 4},
    'sampling_baseline': {'zlib': True, 'complevel': 4},
    'sampling_recent': {'zlib': True, 'complevel': 4},
    'sampling_total': {'zlib': True, 'complevel': 4},
    'continent': {'zlib': True, 'complevel': 4},
}
ds.to_netcdf(out_path, engine='netcdf4', encoding=encoding)
print(f'\nSaved -> {out_path}')
