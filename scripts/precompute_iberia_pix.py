"""Pre-compute Iberian HEALPix-NESTED cell indices for both substrates.

Run once locally (where healpix-geo works cleanly via conda); the two
.npy outputs are small (~few KB each) and committed to the repo.

The DestinE-side download notebook (05_destine_download.py) loads the
nside=128 file at runtime — that way the DestinE platform env does not
need healpix-geo. The Tier-2 clean step (06) uses the nside=64 file as
the analytical grid (matches the Tier-1 HEALPix fit's 110 cells).

Outputs:
    data/precomputed/iberia_pix_nside64_nested.npy   (110 cells)
    data/precomputed/iberia_pix_nside128_nested.npy  (440 cells = 4x parents)

The nside=128 list is the NESTED children of the nside=64 list:
    children = [(parent << 2) | k for k in range(4)]
This guarantees clean parent->child aggregation in Phase D
(no children orphaned, no parents incomplete).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import healpix_geo

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "data" / "precomputed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Must match notebooks/05_destine_download.py:IBERIA_AREA
W, E, S, N = -10.0, 4.0, 35.0, 44.0


# WGS84 (not "sphere"): per the EOPF-DGGS legacy-converters reference
# (legacy_converters/healpix_converters.py:340 + every settings file
# uses {"ellipsoid": {"name": "wgs84"}}), the geo-aware HEALPix
# conversion uses the WGS84 ellipsoid. The sphere↔WGS84 difference is
# small per cell (~hundreds of metres) but compounds across decadal
# climate-impact work — DOMAIN.md flags this explicitly.
ELLIPSOID = "WGS84"


def iberia_nside64() -> np.ndarray:
    """nside=64 NESTED cells whose centres fall in the Iberia bbox.

    Matches the Tier-1 HEALPix analytical grid (Phase C: ~110 cells).
    """
    depth = 6  # log2(64)
    npix = 12 * 64 * 64
    pix = np.arange(npix, dtype=np.int64)
    lons, lats = healpix_geo.nested.healpix_to_lonlat(pix, depth, ELLIPSOID)
    # healpix-geo returns lon in [0, 360); wrap to [-180, 180]
    lons = np.where(lons > 180.0, lons - 360.0, lons)
    mask = (lons >= W) & (lons <= E) & (lats >= S) & (lats <= N)
    return pix[mask]


def main() -> None:
    iberia_64 = iberia_nside64()
    iberia_128 = np.sort(np.concatenate(
        [(iberia_64 << 2) | k for k in range(4)]
    ))

    assert len(iberia_128) == 4 * len(iberia_64)
    parents_recovered = np.unique(iberia_128 >> 2)
    assert np.array_equal(parents_recovered, iberia_64), (
        "parent recovery failed — children/parent invariant broken"
    )

    np.save(OUT_DIR / "iberia_pix_nside64_nested.npy", iberia_64)
    np.save(OUT_DIR / "iberia_pix_nside128_nested.npy", iberia_128)

    print(f"nside=64  (analytical):           {len(iberia_64):>5} cells")
    print(f"nside=128 (download, 4*parents):  {len(iberia_128):>5} cells")
    print(f"Saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
