from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.parse import quote

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response

from .evidence import score_evidence
from .prediction_service import ModelArtifactError, PredictionError, get_predictor
from .schemas import EvidenceRecord, EvidenceRequest, PredictionRequest, PredictionResponse

app = FastAPI(
    title="PhytoTruth AI",
    version="0.2.0",
    description="Evidence-aware natural-products research API. Not for clinical use.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "ptp1b_snapshot.json"
EVIDENCE_INTELLIGENCE_FILE = Path(__file__).resolve().parents[2] / "experiments" / "evidence_intelligence" / "results.json"
KNOWLEDGE_GRAPH_FILE = Path(os.environ["PHYTOTRUTH_KNOWLEDGE_GRAPH"]).expanduser() if os.getenv("PHYTOTRUTH_KNOWLEDGE_GRAPH") else None
DEFAULT_SCREENING_RESULTS_FILE = Path(__file__).resolve().parents[2] / "artifacts" / "npass_screen_100000" / "ranked_candidates.jsonl"
SCREENING_RESULTS_FILE = Path(os.environ["PHYTOTRUTH_SCREENING_RESULTS"]).expanduser() if os.getenv("PHYTOTRUTH_SCREENING_RESULTS") else DEFAULT_SCREENING_RESULTS_FILE


def load_snapshot():
    return json.loads(DATA_FILE.read_text(encoding="utf-8"))


@app.get("/health")
def health():
    return {
        "status": "ok",
        "service": "phytotruth-api",
        "evidence_mode": "experimental + predicted separated",
        "study_planner": EVIDENCE_INTELLIGENCE_FILE.is_file(),
        "knowledge_graph": bool(KNOWLEDGE_GRAPH_FILE and KNOWLEDGE_GRAPH_FILE.is_file()),
        "screening_results": bool(SCREENING_RESULTS_FILE and SCREENING_RESULTS_FILE.is_file()),
    }


@app.get("/api/evidence-intelligence", tags=["Evidence Intelligence"])
def evidence_intelligence():
    """Return the executed study-aware decision-support artifact.

    The artifact is retrospective research planning output, not an autonomous,
    clinical, efficacy, or safety recommendation.
    """
    if not EVIDENCE_INTELLIGENCE_FILE.is_file():
        raise HTTPException(status_code=503, detail="Evidence-intelligence artifact is unavailable")
    return json.loads(EVIDENCE_INTELLIGENCE_FILE.read_text(encoding="utf-8"))


@app.get("/api/dashboard")
def dashboard():
    snapshot = load_snapshot()
    records: list[EvidenceRecord] = []
    molecule_meta = {}
    for molecule in snapshot["molecules"]:
        molecule_meta[molecule["id"]] = {key: value for key, value in molecule.items() if key != "records"}
        for item in molecule["records"]:
            records.append(EvidenceRecord(molecule_id=molecule["id"], **item))
    scores = score_evidence(records)
    candidates = [{**molecule_meta[score.molecule_id], **score.model_dump()} for score in scores]
    potency_order = sorted(candidates, key=lambda item: item["max_pchembl"], reverse=True)
    for index, item in enumerate(potency_order, start=1):
        next(candidate for candidate in candidates if candidate["id"] == item["id"])["potency_rank"] = index
    for index, item in enumerate(candidates, start=1):
        item["evidence_rank"] = index
        item["rank_change"] = item["potency_rank"] - index
    return {
        "metadata": snapshot["metadata"],
        "summary": {
            "molecules": len(candidates),
            "measurements": len(records),
            "independent_assays": len({r.assay_id for r in records}),
            "flagged_candidates": sum(bool(item["flags"]) for item in candidates),
        },
        "candidates": candidates,
    }


@app.post("/api/evidence/score")
def evidence_score(request: EvidenceRequest):
    return score_evidence(request.records, request.active_threshold_pchembl)


@app.post(
    "/api/predict",
    response_model=PredictionResponse,
    tags=["Prediction"],
    summary="Predict target-specific bioactivity from a SMILES structure",
    responses={
        422: {"description": "The request or molecular SMILES is invalid."},
        503: {"description": "The trained model artifact is unavailable or incompatible."},
        500: {"description": "The loaded model could not complete inference."},
    },
)
def predict_smiles(request: PredictionRequest):
    """Run the saved CHEMBL335 model using its training-time RDKit feature configuration.

    The result is a research estimate for a specific target and assay definition. It is
    not experimental confirmation and is not evidence of drug efficacy or clinical safety.
    """
    try:
        predictor = get_predictor()
    except ModelArtifactError as exc:
        raise HTTPException(
            status_code=503,
            detail={"code": "model_unavailable", "message": str(exc)},
        ) from exc
    try:
        return predictor.predict(request.smiles)
    except ValueError as exc:
        raise HTTPException(
            status_code=422,
            detail={"code": "invalid_smiles", "message": "The supplied SMILES is not a valid molecule."},
        ) from exc
    except PredictionError as exc:
        raise HTTPException(
            status_code=500,
            detail={"code": "prediction_failed", "message": str(exc)},
        ) from exc


@app.get("/api/molecule/svg")
def molecule_svg(smiles: str = Query(min_length=1)):
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("invalid smiles")
        drawer = rdMolDraw2D.MolDraw2DSVG(520, 260)
        options = drawer.drawOptions()
        options.clearBackground = False
        options.backgroundColour = (0.02, 0.035, 0.035, 0.0)
        options.bondLineWidth = 1.8
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
        drawer.FinishDrawing()
        return Response(drawer.GetDrawingText(), media_type="image/svg+xml", headers={"Cache-Control": "public, max-age=86400"})
    except ImportError as exc:
        raise HTTPException(503, "RDKit rendering is unavailable") from exc
    except ValueError as exc:
        raise HTTPException(422, "Invalid SMILES") from exc


@app.get("/api/provenance/{molecule_id}")
def provenance(molecule_id: str):
    snapshot = load_snapshot()
    molecule = next((m for m in snapshot["molecules"] if m["id"] == molecule_id), None)
    if not molecule:
        raise HTTPException(404, "Molecule not found")
    return {
        "molecule_id": molecule_id,
        "source_url": f"https://www.ebi.ac.uk/chembl/explore/compound/{quote(molecule_id)}",
        "records": molecule["records"],
        "retrieved_at": snapshot["metadata"]["retrieved_at"],
        "query": snapshot["metadata"]["query"],
    }


@app.get("/api/knowledge-graph/{node_type}/{external_id}", tags=["Knowledge Graph"])
def knowledge_graph_neighborhood(node_type: str, external_id: str, limit: int = Query(default=200, ge=1, le=1000)):
    """Return an evidence-bearing semantic neighborhood, never an unsupported inferred edge."""
    if KNOWLEDGE_GRAPH_FILE is None or not KNOWLEDGE_GRAPH_FILE.is_file():
        raise HTTPException(status_code=503, detail="Knowledge-graph artifact is unavailable")
    from .knowledge_graph import KnowledgeGraph

    with KnowledgeGraph(KNOWLEDGE_GRAPH_FILE, read_only=True) as graph:
        result = graph.neighborhood(node_type, external_id, limit=limit)
    if result is None:
        raise HTTPException(status_code=404, detail="Knowledge-graph node not found")
    return result


@app.get("/api/screening/candidates", tags=["Screening"])
def screening_candidates(limit: int = Query(default=25, ge=1, le=500)):
    """Return precomputed candidates; this endpoint does not trigger an expensive screening run."""
    if SCREENING_RESULTS_FILE is None or not SCREENING_RESULTS_FILE.is_file():
        raise HTTPException(status_code=503, detail="Screening artifact is unavailable")
    candidates = []
    with SCREENING_RESULTS_FILE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                candidates.append(json.loads(line))
                if len(candidates) >= limit:
                    break
    manifest_path = SCREENING_RESULTS_FILE.parent / "screening_manifest.json"
    metadata = None
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        metadata = {
            key: manifest.get(key)
            for key in (
                "created_at", "screened_compounds", "property_filtered_compounds",
                "prediction_failed_compounds", "ranked_candidates", "top_k",
                "out_of_domain_quota", "in_domain_candidates", "out_of_domain_candidates",
                "source_filter", "max_compounds", "disease", "targets", "filters",
                "score_semantics", "limitations", "specifications_sha256",
                "natural_product_database_sha256",
            )
        }
    return {
        "candidates": candidates,
        "count": len(candidates),
        "metadata": metadata,
        "semantics": "research prioritization only; experimental validation is required",
    }

