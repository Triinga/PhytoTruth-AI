from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import logging
import sqlite3
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from .activity import normalize_chembl_activity
from .chembl import iter_activities
from .ingest import NaturalProductIndex, ingest_natural_products
from .structures import standardize_smiles

logger = logging.getLogger(__name__)


def build_index(
    source_file: Path,
    source: str,
    output: Path,
    dataset_version: str,
    member: str | None = None,
    batch_size: int = 1000,
    max_records: int | None = None,
) -> None:
    """Stream a versioned COCONUT/NPASS release into the SQLite provenance index."""
    result = ingest_natural_products(
        source_file,
        source,
        dataset_version,
        output,
        member=member,
        batch_size=batch_size,
        max_records=max_records,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


async def download_chembl(
    target: str,
    output: Path,
    standard_type: str,
    page_size: int = 1000,
    max_records: int | None = None,
) -> None:
    """Stream all ChEMBL pages to JSONL plus a checksum-bearing manifest."""
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".part")
    retrieved_at = datetime.now(UTC).isoformat()
    digest = hashlib.sha256()
    count = 0
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        async for record in iter_activities(
            target,
            standard_type,
            page_size=page_size,
            max_records=max_records,
        ):
            envelope = {
                "source": "ChEMBL",
                "retrieved_at": retrieved_at,
                "target_chembl_id": target,
                "standard_type_filter": standard_type,
                "record": record,
            }
            line = json.dumps(envelope, separators=(",", ":"), ensure_ascii=False) + "\n"
            handle.write(line)
            digest.update(line.encode("utf-8"))
            count += 1
    temporary.replace(output)
    manifest = {
        "source": "ChEMBL REST API",
        "target_chembl_id": target,
        "standard_type": standard_type,
        "retrieved_at": retrieved_at,
        "records": count,
        "page_size": page_size,
        "max_records": max_records,
        "sha256": digest.hexdigest(),
        "data_file": output.name,
    }
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({**manifest, "output": str(output), "manifest": str(manifest_path)}, indent=2))


def _iter_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON on {path}:{line_number}") from exc
            yield line_number, value


def normalize_chembl_file(source: Path, output: Path, rejects: Path) -> dict[str, Any]:
    """Normalize ChEMBL records using an on-disk uniqueness index."""
    output.parent.mkdir(parents=True, exist_ok=True)
    rejects.parent.mkdir(parents=True, exist_ok=True)
    accepted = duplicate_ids = rejected = censored = inconsistent = missing_smiles = 0
    digest = hashlib.sha256()
    output_part = output.with_suffix(output.suffix + ".part")
    rejects_part = rejects.with_suffix(rejects.suffix + ".part")
    with tempfile.TemporaryDirectory(prefix=".phytotruth-", dir=output.parent) as temp_dir:
        connection = sqlite3.connect(Path(temp_dir) / "seen.sqlite")
        connection.execute("CREATE TABLE seen (activity_id TEXT PRIMARY KEY)")
        connection.execute("CREATE TABLE contexts (context_key TEXT PRIMARY KEY)")
        with output_part.open("w", encoding="utf-8", newline="\n") as target, rejects_part.open("w", encoding="utf-8", newline="\n") as bad:
            for line_number, envelope in _iter_jsonl(source):
                raw = envelope.get("record", envelope)
                result = normalize_chembl_activity(raw)
                if result.record is None:
                    bad.write(json.dumps({"source_line": line_number, "reason": result.rejection_reason, "raw": raw}, ensure_ascii=False) + "\n")
                    rejected += 1
                    continue
                cursor = connection.execute("INSERT OR IGNORE INTO seen VALUES (?)", (result.record["source_record_id"],))
                if not cursor.rowcount:
                    duplicate_ids += 1
                    continue
                connection.execute("INSERT OR IGNORE INTO contexts VALUES (?)", (result.record["assay_context_key"],))
                line = json.dumps(result.record, separators=(",", ":"), ensure_ascii=False) + "\n"
                target.write(line)
                digest.update(line.encode("utf-8"))
                censored += int(result.record["is_censored"])
                inconsistent += int(result.record["pchembl_consistent"] is False)
                missing_smiles += int(not result.record.get("canonical_smiles"))
                accepted += 1
            unique_contexts = connection.execute("SELECT COUNT(*) FROM contexts").fetchone()[0]
            connection.close()
    output_part.replace(output)
    rejects_part.replace(rejects)
    summary = {
        "source": "ChEMBL",
        "input": str(source),
        "output": str(output),
        "accepted": accepted,
        "rejected": rejected,
        "duplicate_activity_ids": duplicate_ids,
        "unique_assay_contexts": unique_contexts,
        "censored_measurements": censored,
        "reported_pchembl_inconsistencies": inconsistent,
        "missing_structures": missing_smiles,
        "sha256": digest.hexdigest(),
    }
    output.with_suffix(output.suffix + ".manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def crossmatch_natural_products(source: Path, database: Path, output: Path, rejects: Path) -> dict[str, Any]:
    """Verify ChEMBL molecules against COCONUT/NPASS by standardized InChIKey."""
    output.parent.mkdir(parents=True, exist_ok=True)
    rejects.parent.mkdir(parents=True, exist_ok=True)
    checked = matched = missing_structure = invalid_structure = 0
    matched_sources: set[str] = set()
    digest = hashlib.sha256()
    output_part = output.with_suffix(output.suffix + ".part")
    rejects_part = rejects.with_suffix(rejects.suffix + ".part")
    with NaturalProductIndex(database) as index:
        index.connection.execute("CREATE TEMP TABLE matched_molecules (molecule_id TEXT PRIMARY KEY)")
        with output_part.open("w", encoding="utf-8", newline="\n") as target, rejects_part.open("w", encoding="utf-8", newline="\n") as bad:
            for line_number, record in _iter_jsonl(source):
                checked += 1
                smiles = record.get("canonical_smiles")
                if not smiles:
                    missing_structure += 1
                    bad.write(json.dumps({"source_line": line_number, "reason": "missing canonical_smiles", "record": record}) + "\n")
                    continue
                try:
                    structure = standardize_smiles(smiles)
                except (ValueError, RuntimeError) as exc:
                    invalid_structure += 1
                    bad.write(json.dumps({"source_line": line_number, "reason": str(exc), "record": record}) + "\n")
                    continue
                provenance = index.provenance_for(structure.standard_inchi_key)
                if not provenance:
                    continue
                verified = {
                    **record,
                    "standardized_structure": {
                        "canonical_smiles": structure.canonical_smiles,
                        "standard_inchi": structure.standard_inchi,
                        "standard_inchi_key": structure.standard_inchi_key,
                        "removed_fragments": structure.removed_fragments,
                    },
                    "natural_product_provenance": provenance,
                }
                line = json.dumps(verified, separators=(",", ":"), ensure_ascii=False) + "\n"
                target.write(line)
                digest.update(line.encode("utf-8"))
                matched += 1
                if record.get("molecule_chembl_id"):
                    index.connection.execute("INSERT OR IGNORE INTO matched_molecules VALUES (?)", (record["molecule_chembl_id"],))
                matched_sources.update(item["source"] for item in provenance)
        matched_molecule_count = index.connection.execute("SELECT COUNT(*) FROM matched_molecules").fetchone()[0]
    output_part.replace(output)
    rejects_part.replace(rejects)
    summary = {
        "checked_measurements": checked,
        "matched_measurements": matched,
        "matched_molecules": matched_molecule_count,
        "matched_source_databases": sorted(matched_sources),
        "missing_structure": missing_structure,
        "invalid_structure": invalid_structure,
        "natural_product_index": str(database),
        "sha256": digest.hexdigest(),
    }
    output.with_suffix(output.suffix + ".manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def train(csv_path: Path, model_path: Path) -> None:
    import joblib
    import numpy as np
    import pandas as pd
    from sklearn.ensemble import RandomForestRegressor

    from .features import morgan_fingerprint
    from .ml import scaffold_cross_validate

    frame = pd.read_csv(csv_path).dropna(subset=["canonical_smiles", "pchembl_value"])
    metrics = scaffold_cross_validate(frame.canonical_smiles.tolist(), frame.pchembl_value.tolist())
    X = np.stack([morgan_fingerprint(value) for value in frame.canonical_smiles])
    y = frame.pchembl_value.to_numpy(float)
    model = RandomForestRegressor(n_estimators=600, min_samples_leaf=2, max_features="sqrt", n_jobs=-1, random_state=42)
    model.fit(X, y)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "fingerprint": {"radius": 2, "n_bits": 2048}}, model_path)
    print(json.dumps({"folds": [metric.__dict__ for metric in metrics], "model": str(model_path)}, indent=2))


def reversal_experiment(csv_path: Path, output: Path, k: int) -> None:
    """Molecule-held-out evidence reversal from chronologically or randomly split replicates.

    Input columns: molecule_id, pchembl_value, split (train|heldout), assay_id,
    document_id, canonical_smiles, organism. The split must be assigned before this
    command, ideally by publication date. A molecule needs evidence in both partitions.
    """
    import pandas as pd

    from .evidence import score_evidence
    from .reversal import compare_rankings
    from .schemas import EvidenceRecord

    frame = pd.read_csv(csv_path)
    required = {"molecule_id", "pchembl_value", "split", "assay_id"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    train_rows = frame[frame.split == "train"]
    heldout = frame[frame.split == "heldout"].groupby("molecule_id").pchembl_value.median()
    records = [EvidenceRecord(
        molecule_id=str(row.molecule_id), activity_id=index, assay_id=str(row.assay_id),
        document_id=str(getattr(row, "document_id", "")) or None,
        pchembl_value=float(row.pchembl_value), standard_value_nm=10 ** (9 - float(row.pchembl_value)),
        organism=getattr(row, "organism", None),
    ) for index, row in enumerate(train_rows.itertuples(index=False))]
    scored = {item.molecule_id: item for item in score_evidence(records)}
    ids = sorted(set(scored) & set(heldout.index.astype(str)))
    metrics = compare_rankings(
        [scored[item].max_pchembl for item in ids],
        [scored[item].evidence_score for item in ids],
        [float(heldout.loc[item]) for item in ids],
        k,
    )
    result = {"n_molecules": len(ids), "k": min(k, len(ids)), "metrics": metrics}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="PhytoTruth reproducible pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    chembl = commands.add_parser("fetch-chembl")
    chembl.add_argument("--target", required=True)
    chembl.add_argument("--standard-type", default="IC50")
    chembl.add_argument("--output", type=Path, required=True)
    chembl.add_argument("--page-size", type=int, default=1000)
    chembl.add_argument("--max-records", type=int, help="Explicit safety cap; omit to fetch every API page")
    natural = commands.add_parser("index-natural-products")
    natural.add_argument("--source", choices=["COCONUT", "NPASS"], required=True)
    natural.add_argument("--input", type=Path, required=True)
    natural.add_argument("--output", type=Path, required=True)
    natural.add_argument("--version", required=True, help="Upstream release/version identifier")
    natural.add_argument("--member", help="Exact member for multi-file ZIP archives")
    natural.add_argument("--batch-size", type=int, default=1000)
    natural.add_argument("--max-records", type=int, help="Optional smoke-test cap; omit for the complete release")
    normalize = commands.add_parser("normalize-chembl")
    normalize.add_argument("--input", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    normalize.add_argument("--rejects", type=Path, required=True)
    crossmatch = commands.add_parser("crossmatch-natural-products")
    crossmatch.add_argument("--input", type=Path, required=True, help="Normalized ChEMBL JSONL")
    crossmatch.add_argument("--database", type=Path, required=True)
    crossmatch.add_argument("--output", type=Path, required=True)
    crossmatch.add_argument("--rejects", type=Path, required=True)
    training = commands.add_parser("train")
    training.add_argument("--input", type=Path, required=True)
    training.add_argument("--output", type=Path, default=Path("../models/activity_rf.joblib"))
    activity_training = commands.add_parser("train-activity", help="Run scaffold-held-out RF/XGBoost bioactivity training")
    activity_training.add_argument("--activities", type=Path, required=True, help="Normalized ChEMBL JSONL")
    activity_training.add_argument("--natural-products", type=Path, required=True, help="Natural-product crossmatch JSONL")
    activity_training.add_argument("--output-dir", type=Path, required=True)
    activity_training.add_argument("--context-key")
    activity_training.add_argument("--folds", type=int, default=5)
    activity_training.add_argument("--bootstrap-iterations", type=int, default=1000)
    activity_training.add_argument("--seed", type=int, default=2026)
    reversal = commands.add_parser("reversal")
    reversal.add_argument("--input", type=Path, required=True)
    reversal.add_argument("--output", type=Path, required=True)
    reversal.add_argument("-k", type=int, default=10)
    evidence_reversal = commands.add_parser("evidence-reversal", help="Run the real-data, leakage-resistant ranking benchmark")
    evidence_reversal.add_argument("--activities", type=Path, required=True, help="Normalized ChEMBL JSONL")
    evidence_reversal.add_argument("--natural-products", type=Path, required=True, help="Natural-product crossmatch JSONL")
    evidence_reversal.add_argument("--output-dir", type=Path, required=True)
    evidence_reversal.add_argument("--context-key")
    evidence_reversal.add_argument("--global-cutoff-year", type=int, default=2015)
    evidence_reversal.add_argument("--k", type=int, default=10)
    evidence_reversal.add_argument("--active-threshold", type=float, default=6.0)
    evidence_reversal.add_argument("--bootstrap-iterations", type=int, default=2000)
    evidence_reversal.add_argument("--permutation-iterations", type=int, default=2000)
    evidence_reversal.add_argument("--seed", type=int, default=2026)
    multi_target = commands.add_parser("multi-target-benchmark", help="Run the leakage-resistant multi-target evidence benchmark")
    multi_target.add_argument("--specs", type=Path, required=True)
    multi_target.add_argument("--natural-product-database", type=Path, required=True)
    multi_target.add_argument("--output-dir", type=Path, required=True)
    multi_target.add_argument("--development-fraction", type=float, default=0.6)
    multi_target.add_argument("--k", type=int, default=10)
    multi_target.add_argument("--active-threshold", type=float, default=6.0)
    multi_target.add_argument("--iterations", type=int, default=1000)
    multi_target.add_argument("--seed", type=int, default=2026)
    disagreement = commands.add_parser(
        "assay-disagreement", help="Evaluate prediction of independent later-study assay disagreement"
    )
    disagreement.add_argument("--specs", type=Path, required=True)
    disagreement.add_argument("--output-dir", type=Path, required=True)
    disagreement.add_argument("--bootstrap-iterations", type=int, default=1000)
    disagreement.add_argument("--prediction-interval-iterations", type=int, default=200)
    disagreement.add_argument("--seed", type=int, default=2026)
    generalization = commands.add_parser(
        "assay-disagreement-generalization",
        help="Run study-grouped LOSO selection and locked assay-disagreement evaluation",
    )
    generalization.add_argument("--specs", type=Path, required=True)
    generalization.add_argument("--output-dir", type=Path, required=True)
    generalization.add_argument("--bootstrap-iterations", type=int, default=1000)
    generalization.add_argument("--prediction-interval-iterations", type=int, default=200)
    generalization.add_argument("--seed", type=int, default=2026)
    assay_metadata = commands.add_parser(
        "fetch-assay-metadata",
        help="Fetch exact ChEMBL assay records referenced by disagreement events",
    )
    assay_metadata.add_argument("--events", type=Path, required=True)
    assay_metadata.add_argument("--output", type=Path, required=True)
    assay_metadata.add_argument("--batch-size", type=int, default=50)
    assay_metadata.add_argument("--max-retries", type=int, default=5)
    protocol = commands.add_parser(
        "enrich-assay-protocols",
        help="Extract protocol metadata and run leakage-labelled descriptive analysis",
    )
    protocol.add_argument("--events", type=Path, required=True)
    protocol.add_argument("--assays", type=Path, required=True)
    protocol.add_argument("--output-dir", type=Path, required=True)
    protocol.add_argument("--bootstrap-iterations", type=int, default=1000)
    protocol.add_argument("--seed", type=int, default=2026)
    independent = commands.add_parser(
        "independent-assay-validation",
        help="Score frozen disagreement models on leakage-excluded unseen targets and studies",
    )
    independent.add_argument("--specs", type=Path, required=True)
    independent.add_argument("--protocol", type=Path, required=True)
    independent.add_argument("--previous-specs", type=Path, required=True)
    independent.add_argument("--previous-events", type=Path, required=True)
    independent.add_argument("--full-model", type=Path, required=True)
    independent.add_argument("--selected-model", type=Path, required=True)
    independent.add_argument("--output-dir", type=Path, required=True)
    independent.add_argument("--bootstrap-iterations", type=int, default=2000)
    independent.add_argument("--seed", type=int, default=2026)
    graph = commands.add_parser(
        "build-knowledge-graph",
        help="Build an evidence-traceable plant/organism-compound-target-disease-assay graph",
    )
    graph.add_argument("--specs", type=Path, required=True)
    graph.add_argument("--natural-product-database", type=Path, required=True)
    graph.add_argument("--output", type=Path, required=True)
    screening = commands.add_parser(
        "screen-natural-products",
        help="Run scalable multi-target screening over the NPASS/COCONUT structure index",
    )
    screening.add_argument("--specs", type=Path, required=True)
    screening.add_argument("--natural-product-database", type=Path, required=True)
    screening.add_argument("--output-dir", type=Path, required=True)
    screening.add_argument("--knowledge-graph", type=Path)
    screening.add_argument("--source", choices=["NPASS", "COCONUT"])
    screening.add_argument("--max-compounds", type=int, help="Optional explicit smoke-test cap")
    screening.add_argument("--top-k", type=int, default=100)
    screening.add_argument(
        "--out-of-domain-quota", type=int, default=0,
        help="Reserve this many shortlist positions for the best out-of-domain candidates",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
    if args.command in {"index-natural-products", "crossmatch-natural-products"}:
        from rdkit import RDLogger

        # Structure outcomes are captured as provenance/reject records; routine
        # standardization warnings otherwise overwhelm large-release logs.
        RDLogger.DisableLog("rdApp.warning")
    if args.command == "fetch-chembl":
        asyncio.run(download_chembl(args.target, args.output, args.standard_type, args.page_size, args.max_records))
    elif args.command == "index-natural-products":
        build_index(
            args.input,
            args.source,
            args.output,
            args.version,
            args.member,
            args.batch_size,
            args.max_records,
        )
    elif args.command == "normalize-chembl":
        print(json.dumps(normalize_chembl_file(args.input, args.output, args.rejects), indent=2))
    elif args.command == "crossmatch-natural-products":
        print(json.dumps(crossmatch_natural_products(args.input, args.database, args.output, args.rejects), indent=2))
    elif args.command == "train": train(args.input, args.output)
    elif args.command == "train-activity":
        from .ml_experiment import run_training_experiment

        result = run_training_experiment(
            args.activities,
            args.natural_products,
            args.output_dir,
            context_key=args.context_key,
            folds=args.folds,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "dataset": result["dataset"],
            "selection": result["selection"],
        }, indent=2))
    elif args.command == "reversal": reversal_experiment(args.input, args.output, args.k)
    elif args.command == "evidence-reversal":
        from .reversal_experiment import run_experiment

        result = run_experiment(
            args.activities,
            args.natural_products,
            args.output_dir,
            context_key=args.context_key,
            global_cutoff_year=args.global_cutoff_year,
            k=args.k,
            active_threshold=args.active_threshold,
            bootstrap_iterations=args.bootstrap_iterations,
            permutation_iterations=args.permutation_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "scientific_conclusion": result["scientific_conclusion"],
            "cohort_sizes": {
                name: cohort["evaluation"]["n_molecules"]
                for name, cohort in result["cohorts"].items()
            },
        }, indent=2))
    elif args.command == "multi-target-benchmark":
        from .multi_target_benchmark import run_benchmark

        result = run_benchmark(
            args.specs,
            args.natural_product_database,
            args.output_dir,
            development_fraction=args.development_fraction,
            k=args.k,
            active_threshold=args.active_threshold,
            iterations=args.iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "included_targets": result["included_targets"],
            "selected_ranking": result["development_selection"]["selected_method"],
        }, indent=2))
    elif args.command == "assay-disagreement":
        from .assay_disagreement import run_experiment as run_disagreement_experiment

        result = run_disagreement_experiment(
            args.specs,
            args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations,
            prediction_interval_iterations=args.prediction_interval_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "dataset": result["dataset"],
            "validation_decision": result["validation_decision"],
        }, indent=2))
    elif args.command == "assay-disagreement-generalization":
        from .assay_disagreement_generalization import run_generalization_experiment

        result = run_generalization_experiment(
            args.specs,
            args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations,
            prediction_interval_iterations=args.prediction_interval_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "selected_model": result["development_selection"]["selected_model"],
            "validation_decision": result["validation_decision"],
        }, indent=2))
    elif args.command == "fetch-assay-metadata":
        from .assay_protocol import fetch_assay_metadata

        result = asyncio.run(fetch_assay_metadata(
            args.events,
            args.output,
            batch_size=args.batch_size,
            max_retries=args.max_retries,
        ))
        print(json.dumps(result, indent=2))
    elif args.command == "enrich-assay-protocols":
        from .assay_protocol import run_enrichment

        result = run_enrichment(
            args.events,
            args.assays,
            args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "assays": result["coverage"]["assays"],
            "events": result["coverage"]["events"],
        }, indent=2))
    elif args.command == "independent-assay-validation":
        from .independent_validation import run_independent_validation

        result = run_independent_validation(
            args.specs,
            args.protocol,
            args.previous_specs,
            args.previous_events,
            args.full_model,
            args.selected_model,
            args.output_dir,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        print(json.dumps({
            "output_dir": str(args.output_dir),
            "dataset": result["dataset"],
            "validation_decision": result["validation_decision"],
        }, indent=2))
    elif args.command == "build-knowledge-graph":
        from .knowledge_graph import build_knowledge_graph

        result = build_knowledge_graph(args.specs, args.natural_product_database, args.output)
        print(json.dumps(result, indent=2))
    elif args.command == "screen-natural-products":
        from .screening import screen_natural_products

        result = screen_natural_products(
            args.specs,
            args.natural_product_database,
            args.output_dir,
            knowledge_graph=args.knowledge_graph,
            source=args.source,
            max_compounds=args.max_compounds,
            top_k=args.top_k,
            out_of_domain_quota=args.out_of_domain_quota,
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

