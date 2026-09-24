"""Run only the notebook cells a leave-one-gene-out analysis needs.

Re-executing the whole notebook to reach the LOGO cells at the end costs roughly an
hour, most of it spent on SHAP, figure panels, the deployment forecast, and the
comparison models -- none of which the LOGO analyses depend on. This runner executes
the notebook's *own* cell source, selected by marker string, so the code path is
identical to a full run while taking minutes.

Definition-only cells are truncated at a marker so their expensive bodies are skipped:
the ablation cell contributes ``groups`` / ``drop_groups`` / ``ablation_specs`` without
refitting the 96-model holdout grid, and the first LOGO cell contributes the scoring
helpers without rerunning the leave-one-gene-out sweep.

Nothing here replaces the notebook as the source of truth. Use it while iterating, then
execute the notebook once so its saved outputs carry the final numbers; ``--check``
verifies that the two agree.

Usage:
    python scripts/run_logo_cells.py --out /tmp/logo_fast      # write results elsewhere
    python scripts/run_logo_cells.py --out /tmp/logo_fast --check
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

RELEASE_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = RELEASE_ROOT / "model_essential_nonessential_combined_mutation_level.ipynb"

# (marker identifying the cell, marker at which to truncate it -- None runs the whole cell)
CELL_PLAN: list[tuple[str, str | None]] = [
    ("def find_run_root", None),                       # imports, paths, chdir
    ("catalog_data = pd.read_csv", None),              # load + filter the 2021 features
    ("def create_label_data", None),                   # label_data, numeric_columns, X_*
    ("# === Systematic ablation over feature groups ===", "ablation_results = []"),
    ("# === Gene-grouped (leave-one-gene-out) cross-validation ===", "logo_rows, per_gene_rows = [], []"),
    ("# === Feature-group ablation rescored under gene-grouped", None),
    ("# === Paired bootstrap CIs for single-group ablations", None),
]


def _code_cells() -> list[str]:
    nb = json.loads(NOTEBOOK.read_text())
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def _find(cells: list[str], marker: str) -> str:
    for src in cells:
        if marker in src:
            return src
    raise LookupError(f"no notebook cell contains {marker!r}")


def run(out_dir: Path) -> dict:
    cells = _code_cells()
    ns: dict = {"__name__": "__logo_runner__"}

    # `display` is an IPython builtin the notebook calls; outside IPython it is absent.
    ns["display"] = lambda *a, **k: None

    for marker, cut_at in CELL_PLAN:
        src = _find(cells, marker)
        if cut_at is not None:
            idx = src.find(cut_at)
            if idx == -1:
                raise LookupError(f"truncation marker {cut_at!r} not found in cell {marker!r}")
            src = src[:idx]

        # Redirect the LOGO output directory before the helper cell defines it.
        started = time.time()
        exec(compile(src, f"<notebook:{marker[:40]}>", "exec"), ns)
        if "LOGO_DIR" in ns:
            ns["LOGO_DIR"] = out_dir
            out_dir.mkdir(parents=True, exist_ok=True)
        print(f"  ran {marker[:52]:<52} {time.time() - started:6.1f}s", flush=True)

    return ns


def check(out_dir: Path) -> int:
    """Compare freshly written tables against the notebook's committed ones."""
    import pandas as pd

    committed = RELEASE_ROOT / "results" / "logo"
    failures = 0
    for produced in sorted(out_dir.glob("*.csv")):
        reference = committed / produced.name
        if not reference.exists():
            print(f"  {produced.name}: no committed counterpart yet (new output)")
            continue
        a, b = pd.read_csv(reference), pd.read_csv(produced)
        if a.equals(b):
            print(f"  {produced.name}: matches the notebook")
        else:
            failures += 1
            print(f"  {produced.name}: DIFFERS from the notebook")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True,
                        help="directory for the LOGO tables (keep it outside results/ while iterating)")
    parser.add_argument("--check", action="store_true",
                        help="diff the produced tables against results/logo/")
    args = parser.parse_args()

    out_dir = args.out.resolve()
    started = time.time()
    run(out_dir)
    print(f"\nwrote LOGO tables to {out_dir} in {time.time() - started:.1f}s total")

    if args.check:
        print("\nchecking against the notebook's committed tables:")
        return 1 if check(out_dir) else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
