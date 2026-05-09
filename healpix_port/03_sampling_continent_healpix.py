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

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
HEALPIX_PORT = ROOT / 'healpix_port'
OUT_DIR = HEALPIX_PORT / 'outputs_iberia'
IN_CSV = OUT_DIR / 'bombus_clean_healpix.csv'
PA_NPZ = OUT_DIR / 'presence_absence_healpix.npz'

# ---------------------------------------------------------------------------
# 1. Load cleaned data and the Iberia cell list (the universe of sites).

print('Loading cleaned data + Iberia cell list ...')
df = pd.read_csv(IN_CSV)

pa = np.load(PA_NPZ, allow_pickle=True)
iberia_cells_hp = pa['iberia_cells_hp']
n_cells = int(pa['n_cells'][0])
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

np.savez_compressed(
    OUT_DIR / 'sampling_continent_healpix.npz',
    samp_baseline=samp_baseline,
    samp_recent=samp_recent,
    samp_total=samp_total,
    samp_seasons=samp_stack,
    period_seasons=np.array(period_seasons),
    continent=continent,
    iberia_cells_hp=iberia_cells_hp.astype(np.uint64),
    n_cells=np.array([n_cells], dtype=np.int32),
)
print(f'\nSaved -> {OUT_DIR / "sampling_continent_healpix.npz"}')
