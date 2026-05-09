"""Tier-2 notebooks import this and call ``ensure_destine_or_skip()``.

Detects whether the current process has plausible Destination Earth
Climate Digital Twin (DestinE Climate DT) credentials. If not, the
calling notebook prints a clear notice on stderr and ``sys.exit(0)`` —
so Tier-2 cells in a non-DestinE environment exit cleanly without
breaking the Tier-1 pipeline or the Jupyter Book build.

Detection heuristics (any one match returns True):

  1. ``DESTINE_PLATFORM=true`` env var (manual override / DestinE Jupyter
     image typically sets this).
  2. ``POLYTOPE_USER_KEY`` env var present (matches the polytope-client
     idiom, used by ``earthkit-data`` under the hood for DestinE).
  3. ``~/.polytopeapirc`` exists (the canonical credentials file for
     ``polytope-client``).
  4. ``~/.destine/auth.toml`` exists (Destination Earth Data Lake auth
     bundle written by some DestinE Jupyter images).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def has_destine_credentials() -> bool:
    if os.environ.get("DESTINE_PLATFORM") == "true":
        return True
    if os.environ.get("POLYTOPE_USER_KEY"):
        return True
    if Path("~/.polytopeapirc").expanduser().exists():
        return True
    if Path("~/.destine/auth.toml").expanduser().exists():
        return True
    return False


def ensure_destine_or_skip() -> None:
    if has_destine_credentials():
        return
    print(
        "Tier 2 SKIPPED — DestinE Climate DT credentials not detected.\n"
        "  This notebook only runs on the DestinE Jupyter platform.\n"
        "  Set DESTINE_PLATFORM=true to bypass auto-detection (e.g. for\n"
        "  testing on a CI runner with cached extracts).",
        file=sys.stderr,
    )
    sys.exit(0)


if __name__ == "__main__":
    # Manual smoke test: print detection status and exit cleanly either way.
    if has_destine_credentials():
        print("DestinE credentials detected — Tier 2 would proceed.")
        sys.exit(0)
    ensure_destine_or_skip()
