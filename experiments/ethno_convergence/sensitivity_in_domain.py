"""Prespecified sensitivity analysis: restrict to compounds inside the model
applicability domain.

Listed in `experiments/ethno_convergence_protocol.json` under
`sensitivity_analyses[0]`, written before any outcome was computed.

It re-runs the primary and negative-control comparisons on `score_in_domain`
(the per-plant median over in-domain compounds only) instead of `score`, and
keeps the prespecified `min_compounds_per_plant` threshold, now applied to the
in-domain compound count.

Run from the repository root with the backend virtualenv:

    backend/.venv/Scripts/python.exe experiments/ethno_convergence/sensitivity_in_domain.py

Writes `sensitivity_in_domain.json` beside this file.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))

from app.ethno_convergence import (  # noqa: E402
    family_clustered_bootstrap,
    hodges_lehmann,
    rank_biserial,
)

HERE = Path(__file__).resolve().parent
PROTOCOL = json.loads((ROOT / "experiments" / "ethno_convergence_protocol.json").read_text())
MIN_COMPOUNDS = PROTOCOL["scoring"]["min_compounds_per_plant"]
ITERATIONS = PROTOCOL["uncertainty"]["bootstrap_iterations"]
SEED = PROTOCOL["scoring"]["seed"]
ALPHA = PROTOCOL["primary_outcome"]["alpha"]


def load_rows() -> list[dict]:
    rows: list[dict] = []
    with (HERE / "plant_level_scores.csv").open(newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle):
            if int(record["n_in_domain"] or 0) < MIN_COMPOUNDS:
                continue
            score = float(record["score_in_domain"])
            if math.isnan(score):
                continue
            rows.append(
                {
                    "species": record["species"],
                    "family": record["family"],
                    "score": score,
                    "is_case": record["is_case"] == "True",
                    "is_negative_control": record["is_negative_control"] == "True",
                }
            )
    return rows


def compare(rows: list[dict], *, flag: str) -> dict:
    from scipy.stats import mannwhitneyu

    scoped = [dict(row, is_case=row[flag]) for row in rows]
    case = [row["score"] for row in scoped if row["is_case"]]
    control = [row["score"] for row in scoped if not row["is_case"]]
    low, high = family_clustered_bootstrap(scoped, iterations=ITERATIONS, seed=SEED)
    return {
        "case_plants": len(case),
        "control_plants": len(control),
        "case_median_score": float(sorted(case)[len(case) // 2]) if case else float("nan"),
        "control_median_score": float(sorted(control)[len(control) // 2]) if control else float("nan"),
        "median_shift_hodges_lehmann": hodges_lehmann(case, control),
        "p_value": float(mannwhitneyu(case, control, alternative="two-sided").pvalue),
        "rank_biserial": rank_biserial(case, control),
        "family_clustered_ci_95": [low, high],
    }


def main() -> None:
    rows = load_rows()
    primary = compare(rows, flag="is_case")
    negative = compare(rows, flag="is_negative_control")
    results = {
        "analysis": "sensitivity_in_domain",
        "protocol_version": PROTOCOL["protocol_version"],
        "prespecified_as": "sensitivity_analyses[0]",
        "scope": (
            "Per-plant median over in-domain compounds only; plants retained when they have at "
            f"least {MIN_COMPOUNDS} in-domain compounds."
        ),
        "analysed_plants": len(rows),
        "primary_outcome": primary,
        "negative_control": negative,
        "conditions": {
            "p_below_alpha": primary["p_value"] < ALPHA,
            "ci_excludes_zero": primary["family_clustered_ci_95"][0] > 0
            or primary["family_clustered_ci_95"][1] < 0,
            "negative_control_smaller": abs(negative["median_shift_hodges_lehmann"])
            < abs(primary["median_shift_hodges_lehmann"]),
        },
    }
    results["convergence_supported"] = all(results["conditions"].values())
    (HERE / "sensitivity_in_domain.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
