"""Ethnopharmacological convergence test (protocol: ethno-convergence-v1).

Tests whether plants with a documented traditional antidiabetic use yield
natural products with higher model-predicted PTP1B activity than Duke-covered
control plants.

The activity model never observes ethnobotanical data, so traditional-use
status is an independent label. This measures convergence between two
independent knowledge systems. It establishes no efficacy or safety claim.

All thresholds, the support rule and the negative control are fixed in
``experiments/ethno_convergence_protocol.json`` and were frozen before any
outcome was computed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import random
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import numpy as np

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

NA_TOKENS = {"", "n.a.", "na", "none", "null", "-"}


# --------------------------------------------------------------------------
# normalisation helpers
# --------------------------------------------------------------------------
def normalize_taxon(value: str | None) -> str:
    """Lowercase, collapse whitespace, drop obvious null tokens.

    Deliberately conservative: no synonym resolution and no fuzzy matching, so
    a match is always an exact normalized string identity. This undercounts
    true overlap, which biases the study toward the null.
    """
    if value is None:
        return ""
    cleaned = " ".join(str(value).replace("_", " ").split()).strip().lower()
    return "" if cleaned in NA_TOKENS else cleaned


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matches_any(text: str | None, terms: Sequence[str]) -> bool:
    low = (text or "").lower()
    return any(term in low for term in terms)


# --------------------------------------------------------------------------
# source readers
# --------------------------------------------------------------------------
def _iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


@dataclass(frozen=True)
class DukeUse:
    taxon: str
    activity: str
    family: str
    country: str
    reference: str


def read_duke_ethnobotany(zip_path: Path, member: str = "ETHNOBOT.csv") -> list[DukeUse]:
    """Read Duke ethnobotanical use records. Attribution fields are retained."""
    with zipfile.ZipFile(Path(zip_path)) as archive:
        raw = archive.read(member).decode("utf-8", errors="replace")
    uses: list[DukeUse] = []
    for row in csv.DictReader(io.StringIO(raw)):
        taxon = normalize_taxon(row.get("TAXON"))
        if not taxon:
            continue
        uses.append(
            DukeUse(
                taxon=taxon,
                activity=(row.get("ACTIVITY") or "").strip(),
                family=(row.get("FAMILY") or "").strip(),
                country=(row.get("COUNTRY") or "").strip(),
                reference=(row.get("REFERENCE") or "").strip(),
            )
        )
    return uses


def read_npass_species(path: Path) -> dict[str, dict[str, str]]:
    """org_id -> taxonomy. Falls back to org_name when species_name is absent."""
    organisms: dict[str, dict[str, str]] = {}
    for row in _iter_tsv(path):
        org_id = (row.get("org_id") or "").strip()
        if not org_id:
            continue
        species = normalize_taxon(row.get("species_name")) or normalize_taxon(row.get("org_name"))
        if not species:
            continue
        organisms[org_id] = {
            "species": species,
            "genus": normalize_taxon(row.get("genus_name")),
            "family": normalize_taxon(row.get("family_name")),
            "kingdom": normalize_taxon(row.get("kingdom_name")),
        }
    return organisms


def read_npass_compound_edges(path: Path, keep_orgs: set[str]) -> dict[str, set[str]]:
    """org_id -> {np_id}, streamed; only organisms of interest are retained."""
    edges: dict[str, set[str]] = defaultdict(set)
    for row in _iter_tsv(path):
        org_id = (row.get("org_id") or "").strip()
        if org_id in keep_orgs:
            np_id = (row.get("np_id") or "").strip()
            if np_id:
                edges[org_id].add(np_id)
    return dict(edges)


def read_npass_structures(path: Path, keep_ids: set[str]) -> dict[str, str]:
    """np_id -> SMILES, streamed; only compounds of interest are retained."""
    structures: dict[str, str] = {}
    for row in _iter_tsv(path):
        np_id = (row.get("np_id") or "").strip()
        if np_id not in keep_ids:
            continue
        smiles = ""
        for column in ("canonical_smiles", "SMILES", "smiles", "isomeric_smiles"):
            candidate = (row.get(column) or "").strip()
            if candidate and candidate.lower() not in NA_TOKENS:
                smiles = candidate
                break
        if smiles:
            structures[np_id] = smiles
    return structures


# --------------------------------------------------------------------------
# cohort construction
# --------------------------------------------------------------------------
@dataclass
class Plant:
    species: str
    family: str
    org_ids: set[str] = field(default_factory=set)
    compound_ids: set[str] = field(default_factory=set)
    traditional_uses: list[DukeUse] = field(default_factory=list)
    is_case: bool = False
    is_negative_control: bool = False


def build_plant_cohort(
    *,
    duke_uses: Iterable[DukeUse],
    organisms: dict[str, dict[str, str]],
    diabetes_terms: Sequence[str],
    negative_control_terms: Sequence[str],
) -> dict[str, Plant]:
    """Join NPASS organisms to Duke-covered taxa and assign groups.

    Group assignment uses only traditional-use text: no model output is
    consulted, so the label is independent of the outcome being measured.
    """
    uses_by_taxon: dict[str, list[DukeUse]] = defaultdict(list)
    for use in duke_uses:
        uses_by_taxon[use.taxon].append(use)

    plants: dict[str, Plant] = {}
    for org_id, taxonomy in organisms.items():
        species = taxonomy["species"]
        records = uses_by_taxon.get(species)
        if not records:
            continue  # not Duke-covered: cannot be labelled, so excluded
        plant = plants.get(species)
        if plant is None:
            plant = Plant(
                species=species,
                family=taxonomy.get("family", "") or (records[0].family or "").lower(),
                traditional_uses=records,
                is_case=any(matches_any(r.activity, diabetes_terms) for r in records),
                is_negative_control=any(
                    matches_any(r.activity, negative_control_terms) for r in records
                ),
            )
            plants[species] = plant
        plant.org_ids.add(org_id)
    return plants


def sample_compounds(compound_ids: Sequence[str], cap: int, rng: random.Random) -> list[str]:
    ordered = sorted(compound_ids)  # determinism before sampling
    if len(ordered) <= cap:
        return ordered
    return sorted(rng.sample(ordered, cap))


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def rank_biserial(case: Sequence[float], control: Sequence[float]) -> float:
    """Rank-biserial correlation derived from the Mann-Whitney U statistic."""
    from scipy.stats import mannwhitneyu

    if not case or not control:
        return float("nan")
    u = mannwhitneyu(case, control, alternative="two-sided").statistic
    return float(2.0 * u / (len(case) * len(control)) - 1.0)


def hodges_lehmann(case: Sequence[float], control: Sequence[float]) -> float:
    """Median of all pairwise differences: a robust location shift estimate."""
    if not case or not control:
        return float("nan")
    diffs = np.subtract.outer(np.asarray(case, float), np.asarray(control, float)).ravel()
    return float(np.median(diffs))


def family_clustered_bootstrap(
    rows: Sequence[dict[str, Any]],
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    """Resample whole botanical families, since plants within one are related.

    Returns a 95% percentile interval for the Hodges-Lehmann shift.
    """
    by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_family[row["family"] or f"__unknown__{row['species']}"].append(row)
    families = sorted(by_family)
    if len(families) < 2:
        return (float("nan"), float("nan"))

    rng = np.random.default_rng(seed)
    shifts: list[float] = []
    for _ in range(iterations):
        picked = rng.integers(0, len(families), size=len(families))
        case: list[float] = []
        control: list[float] = []
        for index in picked:
            for row in by_family[families[index]]:
                (case if row["is_case"] else control).append(row["score"])
        if case and control:
            shifts.append(hodges_lehmann(case, control))
    if not shifts:
        return (float("nan"), float("nan"))
    return (float(np.percentile(shifts, 2.5)), float(np.percentile(shifts, 97.5)))


# --------------------------------------------------------------------------
# experiment
# --------------------------------------------------------------------------
def summarise_group(rows: Sequence[dict[str, Any]], *, flag: str) -> dict[str, Any]:
    case = [r["score"] for r in rows if r[flag]]
    control = [r["score"] for r in rows if not r[flag]]
    if not case or not control:
        return {
            "case_plants": len(case),
            "control_plants": len(control),
            "insufficient_data": True,
        }
    from scipy.stats import mannwhitneyu

    test = mannwhitneyu(case, control, alternative="two-sided")
    return {
        "case_plants": len(case),
        "control_plants": len(control),
        "case_median_score": float(np.median(case)),
        "control_median_score": float(np.median(control)),
        "median_shift_hodges_lehmann": hodges_lehmann(case, control),
        "mann_whitney_u": float(test.statistic),
        "p_value": float(test.pvalue),
        "rank_biserial": rank_biserial(case, control),
        "insufficient_data": False,
    }


def run_experiment(
    *,
    species_info: Path,
    species_pair: Path,
    structures: Path,
    duke_zip: Path,
    protocol: dict[str, Any],
    predictor: Any,
    output_dir: Path,
    limit_plants: int | None = None,
) -> dict[str, Any]:
    scoring = protocol["scoring"]
    seed = int(scoring["seed"])
    cap = int(scoring["max_compounds_per_plant"])
    minimum = int(scoring["min_compounds_per_plant"])
    rng = random.Random(seed)

    duke_uses = read_duke_ethnobotany(duke_zip)
    organisms = read_npass_species(species_info)
    plants = build_plant_cohort(
        duke_uses=duke_uses,
        organisms=organisms,
        diabetes_terms=protocol["group_assignment"]["diabetes_activity_terms"],
        negative_control_terms=protocol["negative_control"]["control_activity_terms"],
    )

    org_to_plant = {org: p.species for p in plants.values() for org in p.org_ids}
    edges = read_npass_compound_edges(species_pair, set(org_to_plant))
    for org_id, np_ids in edges.items():
        plants[org_to_plant[org_id]].compound_ids.update(np_ids)

    eligible = sorted(
        (p for p in plants.values() if len(p.compound_ids) >= minimum),
        key=lambda p: p.species,
    )
    if limit_plants is not None:
        eligible = eligible[:limit_plants]

    selection: dict[str, list[str]] = {
        p.species: sample_compounds(sorted(p.compound_ids), cap, rng) for p in eligible
    }
    wanted = {cid for ids in selection.values() for cid in ids}
    smiles_by_id = read_npass_structures(structures, wanted)

    plant_rows: list[dict[str, Any]] = []
    compound_rows: list[dict[str, Any]] = []
    invalid = 0
    for plant in eligible:
        scores: list[float] = []
        in_domain: list[float] = []
        for np_id in selection[plant.species]:
            smiles = smiles_by_id.get(np_id)
            if not smiles:
                invalid += 1
                continue
            try:
                prediction = predictor.predict(smiles)
            except Exception:
                invalid += 1
                continue
            value = float(prediction["predicted_pchembl"])
            inside = bool(prediction["applicability_domain"]["inside"])
            scores.append(value)
            if inside:
                in_domain.append(value)
            compound_rows.append(
                {
                    "species": plant.species,
                    "family": plant.family,
                    "np_id": np_id,
                    "predicted_pchembl": value,
                    "inside_applicability_domain": inside,
                    "is_traditional_antidiabetic": plant.is_case,
                }
            )
        if len(scores) < minimum:
            continue
        diabetes_records = [
            u
            for u in plant.traditional_uses
            if matches_any(u.activity, protocol["group_assignment"]["diabetes_activity_terms"])
        ]
        plant_rows.append(
            {
                "species": plant.species,
                "family": plant.family,
                "is_case": plant.is_case,
                "is_negative_control": plant.is_negative_control,
                "n_compounds_scored": len(scores),
                "score": float(np.median(scores)),
                "score_in_domain": float(np.median(in_domain)) if in_domain else float("nan"),
                "n_in_domain": len(in_domain),
                "traditional_activities": "; ".join(
                    sorted({u.activity for u in diabetes_records})
                ),
                "duke_references": "; ".join(sorted({u.reference for u in diabetes_records if u.reference})),
                "duke_countries": "; ".join(sorted({u.country for u in diabetes_records if u.country})),
            }
        )

    primary = summarise_group(plant_rows, flag="is_case")
    negative = summarise_group(plant_rows, flag="is_negative_control")
    low, high = family_clustered_bootstrap(
        plant_rows,
        iterations=int(protocol["uncertainty"]["bootstrap_iterations"]),
        seed=seed,
    )
    primary["family_clustered_ci_95"] = [low, high]

    conditions = {
        "p_below_alpha": (not primary["insufficient_data"])
        and primary["p_value"] < float(protocol["primary_outcome"]["alpha"]),
        "ci_excludes_zero": bool(
            np.isfinite(low) and np.isfinite(high) and (low > 0 or high < 0)
        ),
        "negative_control_smaller": (
            negative.get("insufficient_data", True)
            or abs(negative.get("median_shift_hodges_lehmann", 0.0))
            < abs(primary.get("median_shift_hodges_lehmann", 0.0))
        ),
    }
    supported = all(conditions.values())

    results = {
        "protocol_version": protocol["protocol_version"],
        "question": protocol["question"],
        "cohort": {
            "duke_use_records": len(duke_uses),
            "npass_organisms": len(organisms),
            "duke_covered_plants": len(plants),
            "eligible_plants": len(eligible),
            "analysed_plants": len(plant_rows),
            "families": len({r["family"] for r in plant_rows}),
            "compounds_scored": len(compound_rows),
            "compounds_unusable": invalid,
        },
        "primary_outcome": primary,
        "negative_control": negative,
        "support_rule": {
            "conditions": conditions,
            "convergence_supported": supported,
            "verdict": "SUPPORTED" if supported else "NOT SUPPORTED",
        },
        "inputs_sha256": {
            "species_info": sha256_file(species_info),
            "species_pair": sha256_file(species_pair),
            "structures": sha256_file(structures),
            "duke": sha256_file(duke_zip),
        },
        "declared_limitations": protocol["declared_limitations"],
    }

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    _write_csv(output_dir / "plant_level_scores.csv", plant_rows)
    _write_csv(output_dir / "compound_level_scores.csv", compound_rows)
    return results


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Ethnopharmacological convergence test")
    parser.add_argument("--species-info", required=True, type=Path)
    parser.add_argument("--species-pair", required=True, type=Path)
    parser.add_argument("--structures", required=True, type=Path)
    parser.add_argument("--duke", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--limit-plants", type=int, default=None)
    args = parser.parse_args(argv)

    from .prediction_service import get_predictor

    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    results = run_experiment(
        species_info=args.species_info,
        species_pair=args.species_pair,
        structures=args.structures,
        duke_zip=args.duke,
        protocol=protocol,
        predictor=get_predictor(),
        output_dir=args.output_dir,
        limit_plants=args.limit_plants,
    )
    print(json.dumps(results["support_rule"], indent=2))
    print(json.dumps(results["cohort"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
