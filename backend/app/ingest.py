from __future__ import annotations

import csv
import gzip
import io
import json
import logging
import re
import sqlite3
import zipfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, BinaryIO, TextIO

from rdkit import Chem

from .structures import standardize_smiles

logger = logging.getLogger(__name__)


def _key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


@dataclass(frozen=True)
class SourceSpec:
    source: str
    id_fields: tuple[str, ...]
    smiles_fields: tuple[str, ...]
    name_fields: tuple[str, ...] = ("name", "common_name", "compound_name")
    organism_fields: tuple[str, ...] = ("organism", "species", "source_organism")
    doi_fields: tuple[str, ...] = ("doi", "reference_doi", "publication_doi")
    collection_fields: tuple[str, ...] = ("collection", "source_collection", "database")


SOURCE_SPECS = {
    "COCONUT": SourceSpec(
        source="COCONUT",
        id_fields=("identifier", "coconut_id", "COCONUT_id", "coconutid", "reference_id"),
        smiles_fields=("canonical_smiles", "canonical smiles", "smiles", "canonicalSmiles"),
    ),
    "NPASS": SourceSpec(
        source="NPASS",
        id_fields=("np_id", "npass_id", "NPASS_ID", "compound_id", "natural_product_id", "identifier"),
        smiles_fields=("canonical_smiles", "canonical smiles", "smiles", "canonicalSmiles"),
    ),
}


@dataclass
class IngestionStats:
    read: int = 0
    accepted_provenance: int = 0
    unique_structures_added: int = 0
    duplicate_structures: int = 0
    duplicate_provenance: int = 0
    invalid_structures: int = 0
    missing_required_fields: int = 0


def _lookup(row: dict[str, Any], aliases: tuple[str, ...]) -> str | None:
    normalized = {_key(str(name)): value for name, value in row.items()}
    for alias in aliases:
        value = normalized.get(_key(alias))
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


@contextmanager
def _open_text(path: Path, member: str | None = None) -> Iterator[TextIO]:
    suffix = path.suffix.lower()
    if suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle
        return
    if suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist() if not name.endswith("/")]
            selected = member or (names[0] if len(names) == 1 else None)
            if selected is None or selected not in names:
                raise ValueError("ZIP contains multiple files; pass --member with an exact member name")
            with archive.open(selected) as binary:
                with io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
                    yield handle
        return
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        yield handle


@contextmanager
def _open_binary(path: Path, member: str | None = None) -> Iterator[BinaryIO]:
    suffix = path.suffix.lower()
    if suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            yield handle
        return
    if suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist() if not name.endswith("/")]
            sdf_names = [name for name in names if name.lower().endswith((".sdf", ".sd"))]
            selected = member or (sdf_names[0] if len(sdf_names) == 1 else None)
            if selected is None or selected not in names:
                raise ValueError("ZIP contains multiple candidate files; pass --member")
            with archive.open(selected) as handle:
                yield handle
        return
    with path.open("rb") as handle:
        yield handle


def _is_sdf(path: Path, member: str | None) -> bool:
    candidate = (member or path.name).lower()
    if candidate.endswith(".gz"):
        candidate = candidate[:-3]
    return candidate.endswith((".sdf", ".sd"))


def _stream_delimited(path: Path, member: str | None) -> Iterator[tuple[int, dict[str, str]]]:
    with _open_text(path, member) as handle:
        sample = handle.read(65536)
        if not sample:
            return
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except csv.Error:
            dialect = csv.excel_tab if "\t" in sample.partition("\n")[0] else csv.excel
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("Delimited source has no header row")
        for line_number, row in enumerate(reader, start=2):
            yield line_number, {str(key): value for key, value in row.items() if key is not None}


def _stream_sdf(path: Path, member: str | None) -> Iterator[tuple[int, dict[str, str]]]:
    with _open_binary(path, member) as handle:
        supplier = Chem.ForwardSDMolSupplier(handle, sanitize=False, removeHs=False, strictParsing=False)
        for record_number, molecule in enumerate(supplier, start=1):
            if molecule is None:
                yield record_number, {"_parse_error": "RDKit could not parse SDF record"}
                continue
            properties = molecule.GetPropsAsDict(includePrivate=False, includeComputed=False)
            row = {str(key): str(value) for key, value in properties.items()}
            if molecule.HasProp("_Name"):
                row.setdefault("identifier", molecule.GetProp("_Name"))
                row.setdefault("name", molecule.GetProp("_Name"))
            try:
                Chem.SanitizeMol(molecule)
                row["canonical_smiles"] = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
            except (ValueError, RuntimeError) as exc:
                row["_parse_error"] = str(exc)
            yield record_number, row


def stream_natural_products(
    path: Path,
    source: str,
    *,
    dataset_version: str = "unspecified",
    member: str | None = None,
    max_records: int | None = None,
) -> Iterator[dict[str, Any]]:
    """Stream official COCONUT/NPASS delimited or SDF exports one record at a time."""
    source_name = source.upper()
    if source_name not in SOURCE_SPECS:
        raise ValueError(f"Unsupported source {source!r}; expected COCONUT or NPASS")
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = SOURCE_SPECS[source_name]
    rows = _stream_sdf(path, member) if _is_sdf(path, member) else _stream_delimited(path, member)
    if max_records is not None and max_records < 1:
        raise ValueError("max_records must be positive when supplied")
    for record_number, (source_line, row) in enumerate(rows, start=1):
        if max_records is not None and record_number > max_records:
            break
        yield {
            "source": source_name,
            "dataset_version": dataset_version,
            "source_file": path.name,
            "source_line": source_line,
            "source_id": _lookup(row, spec.id_fields),
            "smiles": _lookup(row, spec.smiles_fields),
            "name": _lookup(row, spec.name_fields),
            "organism": _lookup(row, spec.organism_fields),
            "doi": _lookup(row, spec.doi_fields),
            "collection": _lookup(row, spec.collection_fields),
            "parse_error": row.get("_parse_error"),
            "raw": row,
        }


class NaturalProductIndex:
    """Disk-backed, constant-memory structure deduplication and provenance index."""

    def __init__(self, path: Path, *, batch_size: int = 1000):
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.batch_size = batch_size
        self.connection = sqlite3.connect(path)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self._pending = 0
        self._create_schema()

    def _create_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS compounds (
                compound_id INTEGER PRIMARY KEY,
                standard_inchi_key TEXT NOT NULL UNIQUE,
                standard_inchi TEXT NOT NULL,
                canonical_smiles TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS provenance (
                provenance_id INTEGER PRIMARY KEY,
                compound_id INTEGER NOT NULL REFERENCES compounds(compound_id),
                source TEXT NOT NULL,
                dataset_version TEXT NOT NULL,
                source_id TEXT NOT NULL,
                source_file TEXT NOT NULL,
                source_line INTEGER NOT NULL,
                original_smiles TEXT NOT NULL,
                name TEXT,
                organism TEXT,
                doi TEXT,
                collection_name TEXT,
                removed_fragments INTEGER NOT NULL,
                raw_json TEXT NOT NULL,
                UNIQUE(source, dataset_version, source_id)
            );
            CREATE TABLE IF NOT EXISTS rejects (
                reject_id INTEGER PRIMARY KEY,
                source TEXT NOT NULL,
                dataset_version TEXT NOT NULL,
                source_file TEXT NOT NULL,
                source_line INTEGER NOT NULL,
                source_id TEXT,
                reason TEXT NOT NULL,
                raw_json TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS idx_provenance_compound ON provenance(compound_id);
            CREATE INDEX IF NOT EXISTS idx_provenance_source ON provenance(source, source_id);
            """
        )
        self.connection.commit()

    def close(self) -> None:
        self.connection.commit()
        self.connection.close()

    def __enter__(self) -> "NaturalProductIndex":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _checkpoint(self) -> None:
        self._pending += 1
        if self._pending >= self.batch_size:
            self.connection.commit()
            self._pending = 0

    def _reject(self, record: dict[str, Any], reason: str) -> None:
        self.connection.execute(
            """INSERT INTO rejects
               (source, dataset_version, source_file, source_line, source_id, reason, raw_json)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (
                record["source"], record["dataset_version"], record["source_file"],
                record["source_line"], record.get("source_id"), reason,
                json.dumps(record.get("raw", {}), sort_keys=True, ensure_ascii=False),
            ),
        )
        self._checkpoint()

    def add(self, record: dict[str, Any], stats: IngestionStats) -> None:
        stats.read += 1
        if record.get("parse_error"):
            stats.invalid_structures += 1
            self._reject(record, f"SDF parse error: {record['parse_error']}")
            return
        if not record.get("source_id") or not record.get("smiles"):
            stats.missing_required_fields += 1
            self._reject(record, "missing source identifier or structure")
            return
        try:
            structure = standardize_smiles(record["smiles"])
        except (ValueError, RuntimeError) as exc:
            stats.invalid_structures += 1
            self._reject(record, f"invalid structure: {exc}")
            return

        cursor = self.connection.execute(
            """INSERT OR IGNORE INTO compounds
               (standard_inchi_key, standard_inchi, canonical_smiles) VALUES (?, ?, ?)""",
            (structure.standard_inchi_key, structure.standard_inchi, structure.canonical_smiles),
        )
        if cursor.rowcount:
            stats.unique_structures_added += 1
        else:
            stats.duplicate_structures += 1
        compound_id = self.connection.execute(
            "SELECT compound_id FROM compounds WHERE standard_inchi_key = ?",
            (structure.standard_inchi_key,),
        ).fetchone()["compound_id"]
        cursor = self.connection.execute(
            """INSERT OR IGNORE INTO provenance
               (compound_id, source, dataset_version, source_id, source_file, source_line,
                original_smiles, name, organism, doi, collection_name, removed_fragments, raw_json)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                compound_id, record["source"], record["dataset_version"], record["source_id"],
                record["source_file"], record["source_line"], structure.original_smiles,
                record.get("name"), record.get("organism"), record.get("doi"),
                record.get("collection"), int(structure.removed_fragments),
                json.dumps(record.get("raw", {}), sort_keys=True, ensure_ascii=False),
            ),
        )
        if cursor.rowcount:
            stats.accepted_provenance += 1
        else:
            stats.duplicate_provenance += 1
        self._checkpoint()

    def provenance_for(self, standard_inchi_key: str) -> list[dict[str, Any]]:
        rows = self.connection.execute(
            """SELECT p.source, p.dataset_version, p.source_id, p.source_file, p.source_line,
                      p.name, p.organism, p.doi, p.collection_name
               FROM provenance p JOIN compounds c USING(compound_id)
               WHERE c.standard_inchi_key = ? ORDER BY p.source, p.source_id""",
            (standard_inchi_key,),
        )
        return [dict(row) for row in rows]

    def summary(self) -> dict[str, int]:
        return {
            "unique_structures": self.connection.execute("SELECT COUNT(*) FROM compounds").fetchone()[0],
            "provenance_records": self.connection.execute("SELECT COUNT(*) FROM provenance").fetchone()[0],
            "rejected_records": self.connection.execute("SELECT COUNT(*) FROM rejects").fetchone()[0],
            "coconut_records": self.connection.execute("SELECT COUNT(*) FROM provenance WHERE source='COCONUT'").fetchone()[0],
            "npass_records": self.connection.execute("SELECT COUNT(*) FROM provenance WHERE source='NPASS'").fetchone()[0],
        }


def ingest_natural_products(
    source_file: Path,
    source: str,
    dataset_version: str,
    database: Path,
    *,
    member: str | None = None,
    batch_size: int = 1000,
    max_records: int | None = None,
) -> dict[str, Any]:
    stats = IngestionStats()
    with NaturalProductIndex(database, batch_size=batch_size) as index:
        for record in stream_natural_products(
            source_file,
            source,
            dataset_version=dataset_version,
            member=member,
            max_records=max_records,
        ):
            index.add(record, stats)
        summary = index.summary()
    result = {
        **asdict(stats),
        **summary,
        "database": str(database),
        "max_records": max_records,
    }
    logger.info("Natural-product ingestion complete", extra=result)
    return result

