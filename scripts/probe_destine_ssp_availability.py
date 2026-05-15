"""Probe DestinE Climate DT for SSP scenario availability.

Run this on the DestinE Jupyter platform (or any host with working DestinE
credentials at ~/.polytopeapirc or POLYTOPE_USER_KEY env var). It probes
every (SSP × model) combination in the cartesian product below by issuing
a minimum-size polytope retrieve — a single day, single timestep, single
variable — and reports which combinations the catalogue accepts.

Output: a structured table on stdout AND a JSON summary at
results/destine_ssp_availability.json.

Expected runtime: ~1-3 minutes (≤ 1 MB per probe × N combinations).
Total download: ~10-30 MB into /tmp (cleaned up automatically).

The probe DOES NOT validate that data is dynamically meaningful — only
that the catalogue accepts the request keys and returns a non-empty GRIB.
A successful probe at one date does not guarantee complete coverage
across the full requested horizon. After identifying populated SSPs,
probe the time-axis edges (e.g. 2049-12-31) separately.
"""
from __future__ import annotations

import json
import os
import tempfile
import time
from datetime import datetime
from pathlib import Path

from polytope.api import Client

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

POLYTOPE_ADDRESS = "https://polytope.lumi.apps.dte.destination-earth.eu"
POLYTOPE_COLLECTION = "destination-earth"

# Cartesian product of candidate scenarios × models. Probe each.
CANDIDATE_SSPS = [
    "SSP1-2.6",
    "SSP2-4.5",
    "SSP3-7.0",   # already validated in the canonical nside=64 chain — sanity check
    "SSP5-8.5",
]
CANDIDATE_MODELS = [
    "IFS-NEMO",
    "IFS-FESOM",
    "ICON",
]

# Reference date in the projection window (well inside any plausible
# scenario archive coverage of the 2020s).
PROBE_DATE = "2030-07-01"

# Mid-century probe to detect scenario-archive truncation (e.g. SSP3-7.0
# was populated through 2039 only as of 2026-05; SSP1-2.6 may have
# different coverage).
PROBE_DATE_LATE = "2049-12-31"

# Cheapest variable to retrieve: total precipitation (param 228,
# accumulated) at the first daily timestep (~50 MB per probe at global
# HEALPix nside=128, vs ~200 MB for 4×/day t2m). Matches the param keys
# used in the canonical chain's 05_destine_download.py.
PROBE_PARAM = "228"
PROBE_TIME = "0000"

OUT_JSON = Path(__file__).resolve().parent.parent / "results" / "destine_ssp_availability.json"

MIN_BYTES = 100_000  # a returned GRIB smaller than this likely indicates empty data


# ---------------------------------------------------------------------------
# Probe logic
# ---------------------------------------------------------------------------

def build_minimal_request(*, experiment: str, model: str, date: str) -> dict:
    """Minimal polytope request — matches the canonical chain's
    05_destine_download.build_request() shape exactly except for the
    1-day date window. The `date` field MUST use the `start/to/end`
    range syntax even for a single-day probe; bare strings raise a
    generic 'Polytope error'.
    """
    return {
        "class": "d1",
        "dataset": "climate-dt",
        "activity": "ScenarioMIP",
        "experiment": experiment,
        "expver": "0001",
        "generation": "1",
        "realization": "1",
        "model": model,
        "resolution": "standard",
        "type": "fc",
        "stream": "clte",
        "levtype": "sfc",
        "param": PROBE_PARAM,
        "date": f"{date}/to/{date}",
        "time": PROBE_TIME,
    }


def probe(client: Client, *, experiment: str, model: str, date: str,
          tmp: Path) -> dict:
    """Issue a minimum-size retrieve. Return a status dict.

    Categorises errors so the summary table is informative without
    requiring the user to read every traceback. Specific error messages
    are not stable polytope contracts — pattern-match cautiously.
    """
    request = build_minimal_request(experiment=experiment, model=model, date=date)
    out_grib = tmp / f"{experiment.replace('.', '_')}__{model}__{date}.grib"
    t0 = time.time()
    try:
        client.retrieve(POLYTOPE_COLLECTION, request, output_file=str(out_grib))
        elapsed = time.time() - t0
        size = out_grib.stat().st_size if out_grib.exists() else 0
        if size < MIN_BYTES:
            return {
                "status": "empty",
                "size_bytes": size,
                "elapsed_s": round(elapsed, 1),
                "note": "request accepted but response shorter than threshold",
            }
        return {
            "status": "ok",
            "size_bytes": size,
            "elapsed_s": round(elapsed, 1),
        }
    except Exception as e:
        elapsed = time.time() - t0
        # Keep the full message (up to 800 chars) — polytope's real
        # diagnostic is often on later lines, so splitlines()[0] discards
        # it. Multiline messages render as escaped \n in the JSON.
        msg_full = str(e).strip() if str(e) else type(e).__name__
        msg_full = msg_full[:800]
        lower = msg_full.lower()
        if "not found" in lower or "no data" in lower or "no values" in lower or "no entries" in lower:
            classification = "catalogue_missing"
        elif "auth" in lower or "credential" in lower or "401" in lower or "403" in lower:
            classification = "auth_error"
        elif "timeout" in lower or "timed out" in lower:
            classification = "timeout"
        else:
            classification = "other_error"
        return {
            "status": classification,
            "elapsed_s": round(elapsed, 1),
            "error": msg_full,
            "error_type": type(e).__name__,
        }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"DestinE Climate DT SSP availability probe")
    print(f"  polytope address:    {POLYTOPE_ADDRESS}")
    print(f"  polytope collection: {POLYTOPE_COLLECTION}")
    print(f"  probe date (mid):    {PROBE_DATE}")
    print(f"  probe date (late):   {PROBE_DATE_LATE}")
    print(f"  probe variable:      param={PROBE_PARAM} time={PROBE_TIME}")
    print(f"  probed at:           {datetime.utcnow().isoformat()}Z")
    print()

    client = Client(address=POLYTOPE_ADDRESS)

    results = {
        "metadata": {
            "address": POLYTOPE_ADDRESS,
            "collection": POLYTOPE_COLLECTION,
            "probe_date_mid": PROBE_DATE,
            "probe_date_late": PROBE_DATE_LATE,
            "probed_at_utc": datetime.utcnow().isoformat() + "Z",
        },
        "combinations": [],
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        for experiment in CANDIDATE_SSPS:
            for model in CANDIDATE_MODELS:
                print(f"--- experiment={experiment!r:>11}  model={model!r:>11} ---")
                # First probe: mid-window
                r_mid = probe(client, experiment=experiment, model=model,
                              date=PROBE_DATE, tmp=tmp)
                if r_mid["status"] == "ok":
                    print(f"  mid  ({PROBE_DATE}): OK  "
                          f"{r_mid['size_bytes']:>12,} B  in {r_mid['elapsed_s']}s")
                    # Only probe late if mid succeeded
                    r_late = probe(client, experiment=experiment, model=model,
                                   date=PROBE_DATE_LATE, tmp=tmp)
                    if r_late["status"] == "ok":
                        print(f"  late ({PROBE_DATE_LATE}): OK  "
                              f"{r_late['size_bytes']:>12,} B  in {r_late['elapsed_s']}s")
                    else:
                        print(f"  late ({PROBE_DATE_LATE}): {r_late['status']}  "
                              f"({r_late.get('error', '')[:80]})")
                else:
                    r_late = {"status": "skipped", "note": "mid-window probe failed"}
                    print(f"  mid  ({PROBE_DATE}): {r_mid['status']}  "
                          f"({r_mid.get('error', '')[:80]})")

                results["combinations"].append({
                    "experiment": experiment,
                    "model": model,
                    "mid_window": r_mid,
                    "late_window": r_late,
                })
                print()

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2))
    print(f"Wrote {OUT_JSON}")

    print()
    print("=" * 80)
    print("  Summary")
    print("=" * 80)
    print(f"  {'experiment':<10}  {'model':<10}  {'mid':<10}  {'late':<10}  notes")
    print("-" * 80)
    for c in results["combinations"]:
        mid = c["mid_window"]["status"]
        late = c["late_window"]["status"]
        note = ""
        if mid == "ok" and late == "ok":
            note = "fully populated (probed both ends)"
        elif mid == "ok" and late != "ok":
            note = f"populated mid-window; late-window missing → truncated archive (cf. SSP3-7.0 pattern)"
        elif mid != "ok":
            note = "not populated for this (experiment, model) combination"
        print(f"  {c['experiment']:<10}  {c['model']:<10}  {mid:<10}  {late:<10}  {note}")
    print("=" * 80)


if __name__ == "__main__":
    main()
