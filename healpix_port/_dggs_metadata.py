"""DGGS Zarr Convention v1 metadata helpers (https://github.com/zarr-conventions/dggs).

Single source of truth for the HEALPix + DGGS attributes attached to every
NetCDF output in the HEALPix-substrate branch of this pipeline.

The convention itself is defined for Zarr; we adapt it to NetCDF by
flattening the nested model into prefix-namespaced attribute keys
(NetCDF group/dataset attributes can't be nested dicts — must be
scalars or 1D arrays). The flattening preserves the convention's
information faithfully so a reader can reconstruct the structured
DGGSZarrConvention / Healpix Pydantic models.

Reference: legacy_converters.core.healpix_conventions
            (https://github.com/EOPF-DGGS/legacy-converters)
"""
from __future__ import annotations

# Per legacy_converters/core/healpix_conventions.py
DGGS_UUID = "7b255807-140c-42ca-97f6-7a1cfecdbc38"
DGGS_SCHEMA_URL = (
    "https://raw.githubusercontent.com/zarr-conventions/dggs/refs/tags/v1/schema.json"
)
DGGS_SPEC_URL = "https://github.com/zarr-conventions/dggs/blob/v1/README.md"
DGGS_DESCRIPTION = "Discrete Global Grid Systems convention for zarr"

# WGS84 ellipsoid constants (per zarr-conventions/dggs v1; also CF-1.10 + EPSG:4326)
WGS84_SEMIMAJOR_AXIS = 6378137.0
WGS84_INVERSE_FLATTENING = 298.257223563

# Canonical names for the spatial dim + coordinate variable (DGGS spec defaults)
SPATIAL_DIMENSION = "cells"
COORDINATE = "cell_ids"


def dggs_attrs(refinement_level: int, indexing_scheme: str = "nested") -> dict:
    """Return a flat attribute dict that encodes the DGGS Zarr Convention v1
    HEALPix grid spec, adapted for NetCDF.

    Parameters
    ----------
    refinement_level
        log2(nside). For nside=64, refinement_level=6.
    indexing_scheme
        One of "nested", "ring", "zuniq". Defaults to "nested" — required
        by DOMAIN.md and used everywhere in this project.
    """
    if indexing_scheme not in ("nested", "ring", "zuniq"):
        raise ValueError(f"indexing_scheme must be nested/ring/zuniq, got {indexing_scheme!r}")
    return {
        # DGGSZarrConvention block (NetCDF attrs accept str/int/float only,
        # so the boolean conformance flag is encoded as int 1)
        "dggs_v1": 1,
        "dggs_uuid": DGGS_UUID,
        "dggs_name": "dggs",
        "dggs_schema_url": DGGS_SCHEMA_URL,
        "dggs_spec_url": DGGS_SPEC_URL,
        "dggs_description": DGGS_DESCRIPTION,
        # Healpix grid spec
        "dggs_grid_name": "healpix",
        "dggs_grid_refinement_level": int(refinement_level),
        "dggs_grid_indexing_scheme": indexing_scheme,
        "dggs_grid_ellipsoid_name": "wgs84",
        "dggs_grid_ellipsoid_semimajor_axis": WGS84_SEMIMAJOR_AXIS,
        "dggs_grid_ellipsoid_inverse_flattening": WGS84_INVERSE_FLATTENING,
        "dggs_grid_spatial_dimension": SPATIAL_DIMENSION,
        "dggs_grid_coordinate": COORDINATE,
        # CF grid_mapping pointer + WKT for the underlying lon/lat datum
        "geospatial_lon_units": "degrees_east",
        "geospatial_lat_units": "degrees_north",
        "geospatial_lonlat_convention": "WGS84, lon in [-180, 180]",
    }


def crs_wgs84_var() -> dict:
    """Return the CF grid_mapping container variable spec for the
    underlying WGS84 lat/lon coordinate reference. Use as one of the
    data_vars in an xarray Dataset; reference it from each measurement
    variable via attrs["grid_mapping"] = "crs_wgs84".
    """
    return {
        "grid_mapping_name": "latitude_longitude",
        "long_name": "WGS84 lat/lon datum (EPSG:4326) for cell-centre coordinates",
        "longitude_of_prime_meridian": 0.0,
        "semi_major_axis": WGS84_SEMIMAJOR_AXIS,
        "inverse_flattening": WGS84_INVERSE_FLATTENING,
        "epsg_code": "EPSG:4326",
        "crs_wkt": (
            'GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
            'ELLIPSOID["WGS 84",6378137,298.257223563]],'
            'CS[ellipsoidal,2],AXIS["latitude",north],'
            'AXIS["longitude",east],ID["EPSG",4326]]'
        ),
        "comment": (
            "Analytical / visualisation projection for downstream EU "
            "biodiversity reporting is ETRS89 / LAEA Europe (EPSG:3035) — "
            "the canonical EEA / Natura 2000 / EUNIS / INSPIRE equal-area "
            "grid. This `crs_wgs84` variable describes the underlying "
            "geographic CRS in which cell-centre coordinates are stored."
        ),
    }


# Convenience constants for the project's specific HEALPix substrate
HEALPIX_NSIDE = 64
HEALPIX_REFINEMENT_LEVEL = 6  # log2(64)
HEALPIX_INDEXING_SCHEME = "nested"
PROJECT_DGGS_ATTRS = dggs_attrs(HEALPIX_REFINEMENT_LEVEL, HEALPIX_INDEXING_SCHEME)
