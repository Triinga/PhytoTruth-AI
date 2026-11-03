import json

from app.multi_target_benchmark import (
    assay_transferability,
    hierarchical_macro_inference,
    load_target_specs,
    select_on_development,
)


def ranking_row(molecule: str, median: float, candidate: float, truth: float):
    return {
        "molecule_id": molecule,
        "median_potency": median,
        "phytotruth_full": candidate,
        "heldout_median_pactivity": truth,
    }


def test_hierarchical_inference_is_reproducible_and_target_aware():
    rows = {
        "T1": [ranking_row("A", 7, 7, 7), ranking_row("B", 6, 6, 6), ranking_row("C", 5, 5, 5)],
        "T2": [ranking_row("D", 5, 7, 7), ranking_row("E", 7, 5, 5), ranking_row("F", 6, 6, 6)],
    }
    first = hierarchical_macro_inference(
        rows, "phytotruth_full", k=2, active_threshold=6, iterations=50, seed=7
    )
    second = hierarchical_macro_inference(
        rows, "phytotruth_full", k=2, active_threshold=6, iterations=50, seed=7
    )
    assert first == second
    assert first["kendall_tau"]["difference"] > 0


def test_development_gate_requires_both_prespecified_metrics():
    passing_metric = {"difference_ci95": [0.01, 0.2], "paired_randomization_p": 0.01, "difference": 0.1}
    failing_metric = {"difference_ci95": [-0.01, 0.2], "paired_randomization_p": 0.2, "difference": 0.1}
    comparisons = {
        "phytotruth_full": {"ndcg_at_k": passing_metric, "kendall_tau": failing_metric},
        "lower_confidence_bound": {"ndcg_at_k": failing_metric, "kendall_tau": passing_metric},
    }
    assert select_on_development(comparisons)["selected_method"] == "median_potency"


def test_transferability_keeps_contexts_separate():
    records = []
    for index in range(10):
        for context, offset, organism in (("A", 0.0, "Homo sapiens"), ("B", 1.0, "Mus musculus")):
            records.append({
                "molecule_chembl_id": f"M{index}",
                "assay_context_key": context,
                "computed_pactivity": 5.0 + index / 10 + offset,
                "assay_context": {
                    "target_organism": organism,
                    "assay_type": "B",
                    "bao_format": "BAO",
                },
            })
    result = assay_transferability(records, {f"M{index}" for index in range(10)})
    assert len(result) == 1
    assert result[0]["same_organism"] is False
    assert result[0]["shared_molecules"] == 10
    assert result[0]["median_absolute_delta"] == 1.0


def test_target_specs_are_unique_and_do_not_contain_outcomes(tmp_path):
    path = tmp_path / "targets.json"
    path.write_text(json.dumps([{
        "target_chembl_id": "CHEMBL1",
        "target_name": "Target",
        "activities": str(tmp_path / "activities.jsonl"),
        "natural_products": str(tmp_path / "natural.jsonl"),
    }]), encoding="utf-8")
    specs = load_target_specs(path)
    assert specs[0]["target_chembl_id"] == "CHEMBL1"
    assert "heldout" not in specs[0]
