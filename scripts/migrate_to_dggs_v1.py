"""One-shot migration script: rename HEALPix dim 'cell' → 'cells' and
coord 'cell' → 'cell_ids' across all NetCDF writers and readers in the
HEALPix-substrate branch, aligning with DGGS Zarr Convention v1.

Run once:
    python scripts/migrate_to_dggs_v1.py

Idempotent: re-running on already-migrated code is a no-op.

Files touched:
  healpix_port/02_presence_absence_healpix.py  (writer)
  healpix_port/03_sampling_continent_healpix.py (writer)
  healpix_port/04_climate_tei_pei_healpix.py    (reader + writer)
  healpix_port/05_regression_healpix.py         (reader)
  notebooks/03h_analysis_healpix.py             (reader)
  notebooks/04h_figures_healpix.py              (reader)
  notebooks/06_destine_clean.py                 (reader + writer)
  notebooks/07_projection.py                    (reader + writer)
  notebooks/08_projection_figures.py            (reader)
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

FILES = [
    "healpix_port/02_presence_absence_healpix.py",
    "healpix_port/03_sampling_continent_healpix.py",
    "healpix_port/04_climate_tei_pei_healpix.py",
    "healpix_port/05_regression_healpix.py",
    "notebooks/03h_analysis_healpix.py",
    "notebooks/04h_figures_healpix.py",
    "notebooks/06_destine_clean.py",
    "notebooks/07_projection.py",
    "notebooks/08_projection_figures.py",
]


def migrate_text(text: str) -> str:
    """Apply DGGS-v1 renames. Returns new text (or unchanged if no-op)."""
    out = text

    # Tuple-dim renames in xarray data_var schemas:
    # ('cell',) and ("cell",) → ('cells',) / ("cells",)
    out = re.sub(r"\(\s*['\"]cell['\"]\s*,\s*\)", "('cells',)", out)
    # ('species', 'cell') and equivalents
    out = re.sub(r"\(\s*['\"]species['\"]\s*,\s*['\"]cell['\"]\s*\)",
                 "('species', 'cells')", out)
    # ('period_season', 'species', 'cell')
    out = re.sub(
        r"\(\s*['\"]period_season['\"]\s*,\s*['\"]species['\"]\s*,\s*['\"]cell['\"]\s*\)",
        "('period_season', 'species', 'cells')", out)

    # Reads: ds['cell'] / ds["cell"] → ds['cell_ids']
    # IMPORTANT: don't touch dict-key uses inside coords={...} blocks where
    # the value is a 1-D array — those become coords={'cell_ids': ('cells', ...)}.
    # We handle those explicitly below.
    # First, the safe reads.
    out = re.sub(r"(\.values|\.attrs|\.sizes\b)\[\s*['\"]cell['\"]\s*\]",
                 r"\1['cell_ids']", out)
    # Generic ds['cell'] when followed by .attrs / .values / .sizes — covered
    # above. But we also need to catch ds['cell'] used in expressions like
    # `pa["cell"].values.astype(...)`. Those are caught by the .values pattern.

    # .sizes['cell'] stayed referring to the dim → 'cells'.
    # Re-do this as a separate substitution because .sizes addresses the DIM
    # not the coord variable.
    out = out.replace('.sizes["cell"]', '.sizes["cells"]')
    out = out.replace(".sizes['cell']", ".sizes['cells']")
    out = out.replace("sizes['cell']", "sizes['cells']")
    out = out.replace('sizes["cell"]', 'sizes["cells"]')

    # coords={"cell": cell_idx} → coords={"cell_ids": ("cells", cell_idx)}
    # This is the writer pattern. Match across newlines.
    coord_pat = re.compile(
        r"""['\"]cell['\"]\s*:\s*([A-Za-z_][\w\.\(\)\,\s\-\+\[\]\d]*)(,?)$""",
        re.MULTILINE,
    )

    def coord_repl(m):
        rhs = m.group(1).rstrip()
        comma = m.group(2)
        return f'"cell_ids": ("cells", {rhs}){comma}'

    # Only replace within coords={...} / .assign_coords(...) blocks. To stay
    # safe, only replace 'cell' as a dict key when the line looks like
    # `'cell': <expr>,` AND its enclosing context is a coords/coord block.
    # Heuristic: the surrounding `coords={` is a few lines above.
    lines = out.splitlines(keepends=True)
    in_coords_block = 0
    for i, line in enumerate(lines):
        if "coords={" in line or "coords = {" in line or ".assign_coords(" in line:
            in_coords_block = 1
        if in_coords_block:
            new_line = re.sub(
                r"['\"]cell['\"]\s*:\s*([^,\n]+)(,?)\s*$",
                lambda m: f'"cell_ids": ("cells", {m.group(1).strip()}){m.group(2)}\n'
                if not line.rstrip().endswith(",")
                else f'"cell_ids": ("cells", {m.group(1).strip()}){m.group(2)}',
                line,
            )
            # Only replace if the line truly looks like a coord assignment
            if new_line != line and ('cell"' in line or "cell'" in line):
                lines[i] = new_line
            if line.rstrip().endswith("}") and "}" in line and "{" not in line:
                in_coords_block = 0
            elif "}" in line and line.count("}") >= line.count("{"):
                in_coords_block = 0
    out = "".join(lines)

    return out


def main():
    for rel in FILES:
        p = ROOT / rel
        if not p.exists():
            print(f"  [skip] {rel} (not found)")
            continue
        before = p.read_text()
        after = migrate_text(before)
        if before == after:
            print(f"  [no-op] {rel}")
        else:
            p.write_text(after)
            n_changed = sum(
                1 for l1, l2 in zip(before.splitlines(), after.splitlines())
                if l1 != l2
            )
            print(f"  [edit]  {rel}  ({n_changed} lines changed)")


if __name__ == "__main__":
    main()
