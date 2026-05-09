"""Probe one DestinE GRIB to extract concrete schema info for Phase D.

Run once after `mamba env update -f environment.yml` brings in eccodes.

Usage:
    python scripts/probe_destine_grib.py
    # or to probe a specific file:
    python scripts/probe_destine_grib.py data/destine/raw/destine_2020_2029_t2m.grib
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import eccodes

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GRIB = ROOT / "data" / "destine" / "raw" / "destine_2020_2029_tp.grib"

# Keys we want to inspect on the first 3 messages. Each tuple is
# (key, fallback) so we don't crash if a particular DestinE GRIB
# doesn't expose a given key.
PROBE_KEYS = [
    ("shortName", "<missing>"),
    ("paramId", -1),
    ("name", "<missing>"),
    ("units", "<missing>"),
    ("dataDate", -1),
    ("dataTime", -1),
    ("step", -1),
    ("stepType", "<missing>"),
    ("numberOfDataPoints", -1),
    ("Nside", -1),
    ("orderingConvention", "<missing>"),
    ("gridType", "<missing>"),
    ("typeOfLevel", "<missing>"),
    ("level", -1),
]


def safe_get(gid: int, key: str, fallback):
    try:
        return eccodes.codes_get(gid, key)
    except Exception:
        return fallback


def probe(grib_path: Path) -> None:
    if not grib_path.exists():
        sys.exit(f"GRIB not found: {grib_path}")

    print(f"=== {grib_path}  ({grib_path.stat().st_size:,} bytes) ===\n")

    msg_count = 0
    sample_count = 3
    first_dates_times = []
    last_dates_times = []
    last_3 = []

    with open(grib_path, "rb") as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                if msg_count < sample_count:
                    record = {k: safe_get(gid, k, fb) for k, fb in PROBE_KEYS}
                    vals = eccodes.codes_get_array(gid, "values")
                    record["values_shape"] = vals.shape
                    record["values_dtype"] = str(vals.dtype)
                    record["values_min"] = float(np.min(vals))
                    record["values_max"] = float(np.max(vals))
                    record["values_mean"] = float(np.mean(vals))
                    print(f"[message {msg_count}]")
                    for k, v in record.items():
                        print(f"  {k}: {v}")
                    print()

                if msg_count < 5:
                    first_dates_times.append(
                        (safe_get(gid, "dataDate", -1),
                         safe_get(gid, "dataTime", -1),
                         safe_get(gid, "step", -1))
                    )

                # Always remember the last few — overwrite a 3-deep buffer
                last_3.append(
                    (safe_get(gid, "dataDate", -1),
                     safe_get(gid, "dataTime", -1),
                     safe_get(gid, "step", -1))
                )
                if len(last_3) > 3:
                    last_3.pop(0)
            finally:
                eccodes.codes_release(gid)
            msg_count += 1

    print(f"Total messages: {msg_count}\n")
    print("First 5 (dataDate, dataTime, step):")
    for d in first_dates_times:
        print(f"  {d}")
    print("\nLast 3 (dataDate, dataTime, step):")
    for d in last_3:
        print(f"  {d}")


if __name__ == "__main__":
    grib = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_GRIB
    probe(grib)
