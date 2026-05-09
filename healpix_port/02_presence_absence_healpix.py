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

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
HEALPIX_PORT = ROOT / 'healpix_port'
OUT_DIR = HEALPIX_PORT / 'outputs_iberia'
IN_CSV = OUT_DIR / 'bombus_clean_healpix.csv'

# HEALPix substrate constants (mirror script 01).
NSIDE = 64
DEPTH = 6
NPIX = 12 * NSIDE * NSIDE          # 49,152

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
    """Pix -> (lon, lat) cell-center, lon wrapped to [-180, 180]."""
    nested = _import_healpix_geo_nested()
    ipix = np.asarray(ipix, dtype='uint64')
    lon, lat = nested.healpix_to_lonlat(ipix, DEPTH)
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
# 6. Save.

np.savez_compressed(
    OUT_DIR / 'presence_absence_healpix.npz',
    species=np.array(species_list),
    period_seasons=np.array(period_seasons),
    iberia_cells_hp=iberia_cells_hp.astype(np.uint64),
    iberia_lon=iberia_lon.astype(np.float32),
    iberia_lat=iberia_lat.astype(np.float32),
    pre_baseline_seasons=np.stack([pre['0_1'], pre['0_2'], pre['0_3']]),
    pre_recent_seasons=np.stack([pre['3_1'], pre['3_2'], pre['3_3']]),
    prab_baseline_seasons=np.stack([prab['0_1'], prab['0_2'], prab['0_3']]),
    prab_recent_seasons=np.stack([prab['3_1'], prab['3_2'], prab['3_3']]),
    prab_baseline=prab_baseline,
    prab_recent=prab_recent,
    sprich_baseline=sprich_baseline,
    sprich_recent=sprich_recent,
    nside=np.array([NSIDE], dtype=np.int32),
    depth=np.array([DEPTH], dtype=np.int32),
    n_cells=np.array([n_cells], dtype=np.int32),
)
print(f'\nSaved -> {OUT_DIR / "presence_absence_healpix.npz"}')

# ---------------------------------------------------------------------------
# 7. Quick per-species summary.

print('\nPer-species cell counts (sample of 10 species):')
for i, spp in enumerate(species_list[:10]):
    n_bs = int((prab_baseline[i] == 1).sum())
    n_rc = int((prab_recent[i] == 1).sum())
    print(f'  {spp:<30} baseline={n_bs:4d}  recent={n_rc:4d}')
print(f'  ... ({len(species_list)} species total)')
