import csv
import sqlite3

import pytest

from app.ingest import NaturalProductIndex, ingest_natural_products, stream_natural_products


def test_coconut_and_npass_streaming_deduplicate_structure_and_keep_provenance(tmp_path):
    coconut = tmp_path / "coconut.csv"
    coconut.write_text(
        "identifier,canonical_smiles,name,organism,doi\n"
        "CNP1,CC(=O)O,Acetic acid,Arabidopsis thaliana,10.0000/example\n"
        "CNP_BAD,not-a-smiles,Bad structure,,\n",
        encoding="utf-8",
    )
    npass = tmp_path / "npass.tsv"
    npass.write_text(
        "NPASS_ID\tCanonical SMILES\tCommon Name\tSpecies\n"
        "NPC1\tCC(=O)[O-].[Na+]\tSodium acetate\tEscherichia coli\n",
        encoding="utf-8",
    )
    database = tmp_path / "natural_products.sqlite"

    coconut_stats = ingest_natural_products(coconut, "COCONUT", "test-release", database, batch_size=1)
    npass_stats = ingest_natural_products(npass, "NPASS", "test-release", database, batch_size=1)

    assert coconut_stats["invalid_structures"] == 1
    assert npass_stats["duplicate_structures"] == 1
    assert npass_stats["unique_structures"] == 1
    assert npass_stats["provenance_records"] == 2
    with NaturalProductIndex(database) as index:
        key = index.connection.execute("SELECT standard_inchi_key FROM compounds").fetchone()[0]
        provenance = index.provenance_for(key)
    assert {item["source"] for item in provenance} == {"COCONUT", "NPASS"}
    assert {item["source_id"] for item in provenance} == {"CNP1", "NPC1"}


def test_duplicate_source_record_does_not_duplicate_provenance(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("identifier,canonical_smiles\nCNP1,CCO\nCNP1,CCO\n", encoding="utf-8")
    database = tmp_path / "index.sqlite"
    stats = ingest_natural_products(source, "COCONUT", "v1", database)
    assert stats["accepted_provenance"] == 1
    assert stats["duplicate_provenance"] == 1
    assert stats["provenance_records"] == 1


def test_zip_with_multiple_members_requires_explicit_member(tmp_path):
    import zipfile

    archive = tmp_path / "release.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("one.csv", "identifier,canonical_smiles\nA,CCO\n")
        handle.writestr("two.csv", "identifier,canonical_smiles\nB,CCN\n")
    records = list(stream_natural_products(archive, "COCONUT", dataset_version="v1", member="two.csv"))
    assert records[0]["source_id"] == "B"
    assert records[0]["dataset_version"] == "v1"


def test_database_constraints_preserve_source_identity(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("identifier,canonical_smiles\nA,CCO\n", encoding="utf-8")
    database = tmp_path / "index.sqlite"
    ingest_natural_products(source, "COCONUT", "v1", database)
    with sqlite3.connect(database) as connection:
        row = connection.execute("SELECT source, dataset_version, source_id, source_file, source_line FROM provenance").fetchone()
    assert row == ("COCONUT", "v1", "A", "source.csv", 2)


def test_streaming_cap_is_explicit_and_validated(tmp_path):
    source = tmp_path / "coconut.tsv"
    source.write_text("COCONUT_ID\tSMILES\nCNP1\tCCO\nCNP2\tCCN\n", encoding="utf-8")

    records = list(stream_natural_products(source, "COCONUT", max_records=1))
    assert [record["source_id"] for record in records] == ["CNP1"]
    with pytest.raises(ValueError, match="positive"):
        list(stream_natural_products(source, "COCONUT", max_records=0))

