"""``farm`` command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .features import FEATURE_COLUMNS, FEATURE_GROUPS, MissingFeaturesError
from .model import available_models, load_model
from .predict import predict


def _cmd_predict(args: argparse.Namespace) -> int:
    try:
        scored = predict(args.input, model=args.model, threshold=args.threshold)
    except MissingFeaturesError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        scored.to_csv(args.output, index=False)
        n_res = int(scored["predicted_resistant"].sum())
        print(
            f"scored {len(scored)} variants with the {args.model} model; "
            f"{n_res} predicted resistance-associated at threshold "
            f"{scored['threshold'].iloc[0]:.3f} -> {args.output}"
        )
    else:
        scored.to_csv(sys.stdout, index=False)
    return 0


def _cmd_info(args: argparse.Namespace) -> int:
    for name in available_models():
        m = load_model(name)
        h = m.metadata.get("holdout", {})
        print(f"{name}")
        print(f"  family              {m.metadata.get('family')}")
        print(f"  decision threshold  {m.threshold:.3f}")
        print(f"  labeled mutations   {m.metadata.get('n_labeled_mutations')}")
        print(f"  holdout AUC         {h.get('auc')}  (n_test={h.get('n_test')})")
        print(f"  holdout sens/spec   {h.get('sensitivity')} / {h.get('specificity')}")
        print()
    print(f"{len(FEATURE_COLUMNS)} features:")
    for group, cols in FEATURE_GROUPS.items():
        print(f"  {group:<24} {len(cols):>2}  {', '.join(cols)}")
    return 0


def _cmd_features(args: argparse.Namespace) -> int:
    if args.json:
        print(json.dumps({"feature_columns": list(FEATURE_COLUMNS),
                          "feature_groups": {k: list(v) for k, v in FEATURE_GROUPS.items()}},
                         indent=2))
    else:
        print("\n".join(FEATURE_COLUMNS))
    return 0


def _cmd_template(args: argparse.Namespace) -> int:
    cols = ["gene", "mutation_oneletter", *FEATURE_COLUMNS]
    frame = pd.DataFrame(columns=cols)
    if args.output:
        frame.to_csv(args.output, index=False)
        print(f"wrote an empty input template with {len(cols)} columns -> {args.output}")
    else:
        frame.to_csv(sys.stdout, index=False)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="farm",
        description="Forecast antibiotic resistance association for M. tuberculosis variants.",
    )
    parser.add_argument("--version", action="version", version=f"farm {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("predict", help="score a CSV of variants carrying FARM features")
    p.add_argument("input", help="CSV with the 25 FARM feature columns")
    p.add_argument("-o", "--output", help="write results here (default: stdout)")
    p.add_argument("-m", "--model", default="combined", choices=available_models(),
                   help="which model to use (default: combined)")
    p.add_argument("-t", "--threshold", type=float, default=None,
                   help="override the calibrated decision threshold")
    p.set_defaults(func=_cmd_predict)

    p = sub.add_parser("info", help="show the bundled models and their calibration")
    p.set_defaults(func=_cmd_info)

    p = sub.add_parser("features", help="list the required feature columns")
    p.add_argument("--json", action="store_true", help="emit JSON instead of one name per line")
    p.set_defaults(func=_cmd_features)

    p = sub.add_parser("template", help="write an empty input CSV with the right columns")
    p.add_argument("-o", "--output", help="write here (default: stdout)")
    p.set_defaults(func=_cmd_template)

    args = parser.parse_args(argv)
    return args.func(args)
