import json

from app.ingest import ingest_natural_products
from app.pipeline import crossmatch_natural_products, normalize_chembl_file


def test_normalize_then_crossmatch_keeps_measurement_and_np_provenance(tmp_path):
    coconut = tmp_path / "coconut.csv"
    coconut.write_text("identifier,canonical_smiles,name\nCNP_ETHANOL,CCO,Ethanol\n", encoding="utf-8")
    database = tmp_path / "np.sqlite"
    ingest_natural_products(coconut, "COCONUT", "fixture-v1", database)

    raw = tmp_path / "chembl.raw.jsonl"
    raw_record = {
        "activity_id": 123,
        "molecule_chembl_id": "CHEMBL_FIXTURE",
        "canonical_smiles": "OCC",
        "assay_chembl_id": "CHEMBL_A",
        "document_chembl_id": "CHEMBL_D",
        "target_chembl_id": "CHEMBL_T",
        "target_organism": "Homo sapiens",
        "standard_type": "IC50",
        "standard_relation": "=",
        "standard_value": "1000",
        "standard_units": "nM",
        "pchembl_value": "6",
        "assay_type": "B",
        "bao_format": "BAO_TEST",
    }
    raw.write_text(json.dumps({"source": "ChEMBL", "record": raw_record}) + "\n", encoding="utf-8")
    normalized = tmp_path / "normalized.jsonl"
    normalization_rejects = tmp_path / "normalization-rejects.jsonl"
    assert normalize_chembl_file(raw, normalized, normalization_rejects)["accepted"] == 1

    matched = tmp_path / "matched.jsonl"
    match_rejects = tmp_path / "match-rejects.jsonl"
    stats = crossmatch_natural_products(normalized, database, matched, match_rejects)
    assert stats["matched_measurements"] == 1
    result = json.loads(matched.read_text(encoding="utf-8"))
    assert result["activity_id"] == 123
    assert result["natural_product_provenance"][0]["source_id"] == "CNP_ETHANOL"
    assert result["standardized_structure"]["standard_inchi_key"]

