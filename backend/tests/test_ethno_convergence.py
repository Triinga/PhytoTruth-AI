"""Tests for the ethnopharmacological convergence experiment.

Fixtures are tiny and deterministic. No external service is contacted and no
measurement is fabricated: the predictor is a stub whose only job is to make
the statistical and bookkeeping logic observable.
"""
from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path

import pytest

from app import ethno_convergence as ec


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture
def duke_zip(tmp_path: Path) -> Path:
    header = [
        "ETHNO", "ACTIVITY", "GENUS", "SPECIES", "SPAUT", "SPRANK", "SPXNAM",
        "SPXAUT", "FAMILY", "CNAME", "COUNTRY", "REFERENCE", "LONGREF",
        "EFFECTIVE", "TAXON", "TAXAUTHOR", "USERID", "CREATED", "MODIFIED",
    ]
    records = [
        ("Diabetes", "Morus", "alba", "Moraceae", "China", "Duke1992", "Morus alba"),
        ("Hyperglycemia", "Panax", "ginseng", "Araliaceae", "Korea", "Duke1992", "Panax ginseng"),
        ("Toothache", "Syzygium", "aromaticum", "Myrtaceae", "India", "Uphof", "Syzygium aromaticum"),
        ("Fever", "Salvia", "officinalis", "Lamiaceae", "Greece", "Uphof", "Salvia officinalis"),
        ("Cough", "Thymus", "vulgaris", "Lamiaceae", "Spain", "Uphof", "Thymus vulgaris"),
    ]
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(header)
    for index, (activity, genus, species, family, country, ref, taxon) in enumerate(records, 1):
        row = [""] * len(header)
        row[0] = str(index)
        row[1] = activity
        row[2] = genus
        row[3] = species
        row[8] = family
        row[10] = country
        row[11] = ref
        row[14] = taxon
        writer.writerow(row)
    path = tmp_path / "duke.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("ETHNOBOT.csv", buffer.getvalue())
    return path


@pytest.fixture
def species_info(tmp_path: Path) -> Path:
    fields = ["org_id", "org_name", "species_name", "genus_name", "family_name", "kingdom_name"]
    rows = [
        {"org_id": "NPO1", "org_name": "Morus alba", "species_name": "Morus alba",
         "genus_name": "Morus", "family_name": "Moraceae", "kingdom_name": "Viridiplantae"},
        {"org_id": "NPO2", "org_name": "Panax ginseng", "species_name": "Panax ginseng",
         "genus_name": "Panax", "family_name": "Araliaceae", "kingdom_name": "Viridiplantae"},
        {"org_id": "NPO3", "org_name": "Salvia officinalis", "species_name": "Salvia officinalis",
         "genus_name": "Salvia", "family_name": "Lamiaceae", "kingdom_name": "Viridiplantae"},
        {"org_id": "NPO4", "org_name": "Thymus vulgaris", "species_name": "Thymus vulgaris",
         "genus_name": "Thymus", "family_name": "Lamiaceae", "kingdom_name": "Viridiplantae"},
        # present in NPASS but absent from Duke: must be excluded from the cohort
        {"org_id": "NPO5", "org_name": "Obscura nondukea", "species_name": "Obscura nondukea",
         "genus_name": "Obscura", "family_name": "Unknownaceae", "kingdom_name": "Viridiplantae"},
    ]
    return _write_tsv(tmp_path / "species_info.txt", fields, rows)


@pytest.fixture
def species_pair(tmp_path: Path) -> Path:
    fields = ["src_org_record_id", "org_id", "np_id"]
    rows = []
    plan = {"NPO1": 4, "NPO2": 4, "NPO3": 4, "NPO4": 4, "NPO5": 4}
    counter = 0
    for org_id, count in plan.items():
        for _ in range(count):
            counter += 1
            rows.append({"src_org_record_id": str(counter), "org_id": org_id,
                         "np_id": f"NPC{counter:04d}"})
    return _write_tsv(tmp_path / "species_pair.txt", fields, rows)


@pytest.fixture
def structures(tmp_path: Path, species_pair: Path) -> Path:
    fields = ["np_id", "canonical_smiles"]
    rows = []
    with species_pair.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({"np_id": row["np_id"], "canonical_smiles": "CCO"})
    return _write_tsv(tmp_path / "structures.txt", fields, rows)


@pytest.fixture
def protocol() -> dict:
    root = Path(__file__).resolve().parents[2]
    return json.loads(
        (root / "experiments" / "ethno_convergence_protocol.json").read_text(encoding="utf-8")
    )


class StubPredictor:
    """Returns a fixed score per species group so the test outcome is known."""

    def __init__(self, case_ids: set[str], high: float = 7.0, low: float = 5.0) -> None:
        self.case_ids = case_ids
        self.high = high
        self.low = low
        self.calls = 0

    def predict(self, smiles: str) -> dict:
        self.calls += 1
        # score encoded by call order is not used; see keyed stub below
        return {
            "predicted_pchembl": self.low,
            "applicability_domain": {"inside": True},
        }


class KeyedPredictor:
    """Scores by compound id prefix, injected through a SMILES suffix comment."""

    def __init__(self, high_ids: set[str]) -> None:
        self.high_ids = high_ids

    def predict(self, smiles: str) -> dict:
        value = 7.5 if smiles in self.high_ids else 5.0
        return {
            "predicted_pchembl": value,
            "applicability_domain": {"inside": True},
        }


# --------------------------------------------------------------------------
# unit tests
# --------------------------------------------------------------------------
def test_normalize_taxon_handles_nulls_and_whitespace():
    assert ec.normalize_taxon("  Morus   ALBA ") == "morus alba"
    assert ec.normalize_taxon("n.a.") == ""
    assert ec.normalize_taxon(None) == ""
    assert ec.normalize_taxon("Panax_ginseng") == "panax ginseng"


def test_group_assignment_uses_only_traditional_use_text(duke_zip, species_info, protocol):
    plants = ec.build_plant_cohort(
        duke_uses=ec.read_duke_ethnobotany(duke_zip),
        organisms=ec.read_npass_species(species_info),
        diabetes_terms=protocol["group_assignment"]["diabetes_activity_terms"],
        negative_control_terms=protocol["negative_control"]["control_activity_terms"],
    )
    assert plants["morus alba"].is_case is True
    assert plants["panax ginseng"].is_case is True       # matched via "hyperglyc"
    assert plants["salvia officinalis"].is_case is False
    assert plants["thymus vulgaris"].is_case is False


def test_plants_absent_from_duke_are_excluded(duke_zip, species_info, protocol):
    """Conditioning on Duke coverage prevents confounding use with coverage."""
    plants = ec.build_plant_cohort(
        duke_uses=ec.read_duke_ethnobotany(duke_zip),
        organisms=ec.read_npass_species(species_info),
        diabetes_terms=protocol["group_assignment"]["diabetes_activity_terms"],
        negative_control_terms=protocol["negative_control"]["control_activity_terms"],
    )
    assert "obscura nondukea" not in plants


def test_traditional_use_attribution_is_retained(duke_zip, species_info, protocol):
    plants = ec.build_plant_cohort(
        duke_uses=ec.read_duke_ethnobotany(duke_zip),
        organisms=ec.read_npass_species(species_info),
        diabetes_terms=protocol["group_assignment"]["diabetes_activity_terms"],
        negative_control_terms=protocol["negative_control"]["control_activity_terms"],
    )
    uses = plants["morus alba"].traditional_uses
    assert uses[0].reference == "Duke1992"
    assert uses[0].country == "China"


def test_compound_sampling_is_deterministic_and_capped():
    import random

    ids = [f"NPC{i:04d}" for i in range(200)]
    first = ec.sample_compounds(ids, 50, random.Random(2026))
    second = ec.sample_compounds(ids, 50, random.Random(2026))
    assert first == second
    assert len(first) == 50
    assert set(first) <= set(ids)


def test_sampling_returns_all_when_below_cap():
    import random

    ids = ["NPC2", "NPC1", "NPC3"]
    assert ec.sample_compounds(ids, 50, random.Random(1)) == ["NPC1", "NPC2", "NPC3"]


def test_hodges_lehmann_detects_positive_shift():
    assert ec.hodges_lehmann([7.0, 7.5, 8.0], [5.0, 5.5]) == pytest.approx(2.25, abs=0.3)


def test_rank_biserial_sign_follows_direction():
    assert ec.rank_biserial([9.0, 9.5], [1.0, 1.5]) > 0
    assert ec.rank_biserial([1.0, 1.5], [9.0, 9.5]) < 0


def test_structures_reader_keeps_only_requested_ids(structures):
    kept = ec.read_npass_structures(structures, {"NPC0001", "NPC0002"})
    assert set(kept) == {"NPC0001", "NPC0002"}


def test_edges_reader_keeps_only_requested_orgs(species_pair):
    edges = ec.read_npass_compound_edges(species_pair, {"NPO1"})
    assert set(edges) == {"NPO1"}
    assert len(edges["NPO1"]) == 4


# --------------------------------------------------------------------------
# end-to-end
# --------------------------------------------------------------------------
def test_experiment_detects_a_planted_signal(
    tmp_path, duke_zip, species_info, species_pair, structures, protocol
):
    """With case plants scored high by construction, the shift must be positive."""
    case_smiles = {"CASE"}
    # rewrite structures so case-plant compounds carry a distinguishable SMILES
    org_of = {}
    with species_pair.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            org_of[row["np_id"]] = row["org_id"]
    rows = [
        {"np_id": np_id, "canonical_smiles": "CASE" if org in {"NPO1", "NPO2"} else "CTRL"}
        for np_id, org in org_of.items()
    ]
    keyed = _write_tsv(tmp_path / "structures_keyed.txt", ["np_id", "canonical_smiles"], rows)

    local = dict(protocol)
    local["scoring"] = {**protocol["scoring"], "min_compounds_per_plant": 3,
                        "max_compounds_per_plant": 50}
    local["uncertainty"] = {**protocol["uncertainty"], "bootstrap_iterations": 200}

    results = ec.run_experiment(
        species_info=species_info,
        species_pair=species_pair,
        structures=keyed,
        duke_zip=duke_zip,
        protocol=local,
        predictor=KeyedPredictor(case_smiles),
        output_dir=tmp_path / "out",
    )
    primary = results["primary_outcome"]
    assert primary["case_plants"] == 2
    assert primary["control_plants"] == 2
    assert primary["median_shift_hodges_lehmann"] == pytest.approx(2.5, abs=0.01)
    assert results["cohort"]["duke_covered_plants"] == 4
    assert (tmp_path / "out" / "results.json").exists()
    assert (tmp_path / "out" / "plant_level_scores.csv").exists()


def test_experiment_reports_not_supported_without_signal(
    tmp_path, duke_zip, species_info, species_pair, structures, protocol
):
    """Flat scores must not pass the support rule."""
    local = dict(protocol)
    local["scoring"] = {**protocol["scoring"], "min_compounds_per_plant": 3}
    local["uncertainty"] = {**protocol["uncertainty"], "bootstrap_iterations": 200}

    results = ec.run_experiment(
        species_info=species_info,
        species_pair=species_pair,
        structures=structures,
        duke_zip=duke_zip,
        protocol=local,
        predictor=KeyedPredictor(set()),  # every compound scores identically
        output_dir=tmp_path / "flat",
    )
    assert results["support_rule"]["convergence_supported"] is False
    assert results["support_rule"]["verdict"] == "NOT SUPPORTED"


def test_input_checksums_are_recorded(
    tmp_path, duke_zip, species_info, species_pair, structures, protocol
):
    local = dict(protocol)
    local["uncertainty"] = {**protocol["uncertainty"], "bootstrap_iterations": 50}
    results = ec.run_experiment(
        species_info=species_info,
        species_pair=species_pair,
        structures=structures,
        duke_zip=duke_zip,
        protocol=local,
        predictor=KeyedPredictor(set()),
        output_dir=tmp_path / "sums",
    )
    for key in ("species_info", "species_pair", "structures", "duke"):
        assert len(results["inputs_sha256"][key]) == 64
