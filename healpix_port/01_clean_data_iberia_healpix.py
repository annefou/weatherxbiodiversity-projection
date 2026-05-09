# HEALPix-NESTED nside=64 port of soroye_port/01_clean_data_iberia.py — substrate-robustness branch.
"""
HEALPix-NESTED nside=64 port of script 01.

Mirrors the upstream Iberia clean (Kerr-2015 species filter, err-obs
removal, period 0/3 selection, season assignment) but assigns each
occurrence to a HEALPix nside=64 NESTED cell using `healpix-geo`
(geo-aware, WGS84-friendly). The on-disk pixelisation is **always
NESTED** — see `DOMAIN.md` § "HEALPix is always NESTED".

At nside=64 the HEALPix cell area is 41,253/(12·64²) ~= 53.7 deg^2 -- per
the equal-area HEALPix tessellation that's ~92x92 km, scale-matched to
Soroye's 100 km cylindrical-equal-area grid. The substrate-robustness
test in this branch checks whether `sc_TEI_delta` is recovered at the
same magnitude on this independent equal-area substrate.

Output: `healpix_port/outputs_iberia/bombus_clean_healpix.csv` — same
columns as the CEA version PLUS `cell_id_hp` (HEALPix nside=64 NESTED
pixel index, np.uint64).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / 'reference' / 'Bumblebee_repo_wbombusdat'
GBIF_CSV = ROOT / 'data' / 'gbif_dl' / '0006204-260423192947929.csv'
KERR_SPP_CSV = REF / '0_data' / 'Kerr_et_al2015_spplist.csv'
ERR_OBS_CSV = REF / '0_data' / 'bombus_err_obs.csv'

OUT_DIR = ROOT / 'healpix_port' / 'outputs_iberia'
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CLEAN = OUT_DIR / 'bombus_clean_healpix.csv'

# HEALPix substrate constants — depth = log2(nside) for healpix-geo.
NSIDE = 64
DEPTH = 6                          # 2**6 == 64
NPIX = 12 * NSIDE * NSIDE          # 49,152 global cells

# Iberia bbox (matches notebooks/05_destine_download.py:IBERIA_AREA).
IBERIA_LON_MIN, IBERIA_LON_MAX = -10.0, 4.0
IBERIA_LAT_MIN, IBERIA_LAT_MAX = 35.0, 44.0


# ---------------------------------------------------------------------------
# HEALPix helpers — defensive import shim against API drift.

def _import_healpix_geo_nested():
    """Return module healpix_geo.nested. Defensive: try a couple of
    likely import paths so this script tolerates minor API shifts in
    `healpix-geo` 0.1.x. We never fall back to `healpy` because it is
    non-geo-aware (cosmology lon convention 0-360, no datum) and
    accumulates small biases over decadal time series — unacceptable
    for biodiversity-precision climate-impact work (`DOMAIN.md`).
    """
    try:
        from healpix_geo import nested
        return nested
    except Exception:
        # Some older versions exposed the same surface at the top level.
        import healpix_geo as nested
        return nested


def lonlat_to_pix(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Vectorised lon/lat -> HEALPix NESTED pixel index at depth=DEPTH.

    `healpix-geo` works on the WGS84 ellipsoid; we pass `ellipsoid='sphere'`
    here for parity with the cell-center round-trip used downstream
    (HEALPix is intrinsically defined on the sphere; the 'WGS84' option
    in healpix-geo applies a tiny correction we don't need at the
    biodiversity scale).
    """
    nested = _import_healpix_geo_nested()
    lon = np.asarray(lon, dtype='float64')
    lat = np.asarray(lat, dtype='float64')
    return nested.lonlat_to_healpix(lon, lat, DEPTH).astype(np.uint64)


def pix_to_lonlat(ipix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised pix -> (lon, lat) cell-center, lon wrapped to
    [-180, 180]. healpix-geo returns lon in [0, 360) by default."""
    nested = _import_healpix_geo_nested()
    ipix = np.asarray(ipix, dtype='uint64')
    lon, lat = nested.healpix_to_lonlat(ipix, DEPTH)
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    return lon, lat


# ---------------------------------------------------------------------------
# 1. Load raw GBIF — same TSV as the CEA pipeline.

print('Loading GBIF Iberia data ...')
gbif = pd.read_csv(GBIF_CSV, sep='\t', low_memory=False)
print(f'  raw: {len(gbif):,} rows')

gbif = gbif.rename(columns={
    'decimalLatitude':  'latitude',
    'decimalLongitude': 'longitude',
    'year':             'year_',
    'gbifID':           'LYID',
    'stateProvince':    'state_prov',
    'countryCode':      'country',
})
gbif['species'] = (
    gbif['species']
    .astype(str)
    .str.replace(r'^Bombus\s+', '', regex=True)
    .str.lower()
)
gbif['year_'] = pd.to_numeric(gbif['year_'], errors='coerce')

# ---------------------------------------------------------------------------
# 2. Clean (mirror of upstream script 01)

drop_species = {
    'x.virginica', 'megachile_sculpturalis', 'fernaldaepsithyrus',
    'subterraneobombus', 'psithyrus', 'laesobombus', 'allopsithyrus',
}

bombus_clean = (
    gbif
    .dropna(subset=['latitude', 'longitude', 'species'])
    .query('year_ >= 1901 and year_ < 2015')
    .loc[lambda d: ~d['species'].isin(drop_species)]
    .copy()
)
print(f'  after basic cleaning: {len(bombus_clean):,} rows')

# ---------------------------------------------------------------------------
# 3. Continent + period assignment

bombus_clean['continent'] = 2   # Iberia is in Europe


def assign_period(year: float) -> float:
    if 1901 <= year < 1975: return 0
    if 1975 <= year < 1987: return 1
    if 1987 <= year < 2000: return 2
    if 2000 <= year < 2015: return 3
    return np.nan


bombus_clean['timeperiod'] = bombus_clean['year_'].apply(assign_period)

# ---------------------------------------------------------------------------
# 4. Species synonyms

synonyms = {
    'ashtoni': 'bohemicus',
    'californicus': 'fervidus',
    'moderatus': 'cryptarum',
    'sonorus': 'pensylvanicus',
    'alboanalis': 'jonellus',
    'wilmattae': 'ephippiatus',
    'volucelloides': 'melaleucus',
    'soroensis': 'soroeensis',
    'sichelii': 'sicheli',
    'wurfleni': 'wurflenii',
    'bannitus': 'muscorum',
    'mocsaryi': 'laesus',
}
bombus_clean['species'] = bombus_clean['species'].replace(synonyms)

# ---------------------------------------------------------------------------
# 5. Filter to Kerr 2015 species list

kerr = pd.read_csv(KERR_SPP_CSV)
kerr_species = set(kerr['species'].astype(str).str.lower().unique())
print(f'  Kerr species total: {len(kerr_species)}')

bombus_clean = bombus_clean[bombus_clean['species'].isin(kerr_species)].copy()
print(f'  after Kerr filter: {len(bombus_clean):,} rows  ({bombus_clean["species"].nunique()} species)')

# ---------------------------------------------------------------------------
# 6. Remove erroneous obs (apply the same table as upstream)

err_obs = pd.read_csv(ERR_OBS_CSV, na_values=[''])

kept_parts = []
for _, row in err_obs.iterrows():
    species = row['species']
    excl_cont = row['exclude_from_cont']
    excl_state_raw = row['exclude_from_state']
    excl_states = (
        [s.strip() for s in str(excl_state_raw).split(',')]
        if pd.notna(excl_state_raw) else []
    )

    subset = bombus_clean[bombus_clean['species'] == species].copy()
    if pd.notna(excl_cont):
        subset = subset[subset['continent'] != excl_cont]
    if excl_states:
        keep = subset['state_prov'].isna() | ~subset['state_prov'].isin(excl_states)
        subset = subset[keep]
    kept_parts.append(subset)

species_in_err = set(err_obs['species'].astype(str).str.lower().unique())
kept_parts.append(bombus_clean[~bombus_clean['species'].isin(species_in_err)])

bombus_clean = pd.concat(kept_parts, ignore_index=True)
print(f'  after err-obs removal: {len(bombus_clean):,} rows')

# ---------------------------------------------------------------------------
# 7. Rename + add season + drop middle periods

bombus_clean = bombus_clean.rename(columns={'year_': 'year', 'timeperiod': 'period'})
bombus_clean = bombus_clean[bombus_clean['period'].isin([0, 3])].copy()


def assign_season(year: float) -> float:
    if 1901 <= year <= 1924: return 1
    if 1925 <= year <= 1949: return 2
    if 1950 <= year <= 1974: return 3
    if 2000 <= year <= 2004: return 1
    if 2005 <= year <= 2009: return 2
    if 2010 <= year <= 2014: return 3
    return np.nan


bombus_clean['season'] = bombus_clean['year'].apply(assign_season)
bombus_clean['period_season'] = (
    bombus_clean['period'].astype(int).astype(str)
    + '_' + bombus_clean['season'].astype(int).astype(str)
)

# ---------------------------------------------------------------------------
# 8. HEALPix cell assignment — NESTED nside=64.

print(f'\nAssigning HEALPix nside={NSIDE} NESTED cell ids '
      f'(global npix={NPIX:,}) ...')
bombus_clean['cell_id_hp'] = lonlat_to_pix(
    bombus_clean['longitude'].values,
    bombus_clean['latitude'].values,
)

cols = ['species', 'latitude', 'longitude', 'year', 'LYID',
        'continent', 'period', 'season', 'period_season', 'cell_id_hp']
bombus_clean = bombus_clean[cols].reset_index(drop=True)

bombus_clean.to_csv(OUT_CLEAN, index=False)
print(f'\nSaved -> {OUT_CLEAN}  ({len(bombus_clean):,} rows)')

print('\nPer-period/season summary:')
print(bombus_clean.groupby('period_season').agg(
    nspp=('species', 'nunique'),
    nLYID=('LYID', 'nunique'),
    n_records=('species', 'size'),
    n_cells_hp=('cell_id_hp', 'nunique'),
).to_string())
print(f'\nSpecies: {bombus_clean["species"].nunique()}')

# ---------------------------------------------------------------------------
# 9. Diagnostic: how many HEALPix cells does Iberia contain at this depth?
#    (computed independently of the GBIF data — gives the user a sanity
#    check that the substrate is dimensioned correctly).

all_pix = np.arange(NPIX, dtype='uint64')
all_lon, all_lat = pix_to_lonlat(all_pix)
iberia_mask = (
    (all_lon >= IBERIA_LON_MIN) & (all_lon <= IBERIA_LON_MAX)
    & (all_lat >= IBERIA_LAT_MIN) & (all_lat <= IBERIA_LAT_MAX)
)
n_iberia_total = int(iberia_mask.sum())

# Cells with at least one occurrence in baseline + recent.
period_counts = (
    bombus_clean.groupby(['period', 'cell_id_hp']).size().reset_index(name='n')
)
n_cells_bs = int(period_counts.loc[period_counts['period'] == 0, 'cell_id_hp'].nunique())
n_cells_rc = int(period_counts.loc[period_counts['period'] == 3, 'cell_id_hp'].nunique())
n_cells_either = int(bombus_clean['cell_id_hp'].nunique())

print(f'\nHEALPix nside={NSIDE} Iberia coverage:')
print(f'  total Iberian cells (bbox lon {IBERIA_LON_MIN}..{IBERIA_LON_MAX}, '
      f'lat {IBERIA_LAT_MIN}..{IBERIA_LAT_MAX}): {n_iberia_total}')
print(f'  cells with baseline (period 0) occurrences:   {n_cells_bs}')
print(f'  cells with recent  (period 3) occurrences:    {n_cells_rc}')
print(f'  cells with any occurrence (baseline OR recent): {n_cells_either}')
