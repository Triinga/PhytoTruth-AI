from __future__ import annotations

import asyncio
import logging
from collections.abc import AsyncIterator
from typing import Any, Sequence
from urllib.parse import urljoin

import httpx

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data/"
MAX_PAGE_SIZE = 1000
RETRYABLE_STATUS_CODES = {408, 425, 429, 500, 502, 503, 504}
logger = logging.getLogger(__name__)


class ChEMBLRequestError(RuntimeError):
    """Raised after a ChEMBL request exhausts its retry budget."""


def _retry_delay(response: httpx.Response | None, attempt: int, maximum: float) -> float:
    if response is not None:
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                return min(float(retry_after), maximum)
            except ValueError:
                pass
    return min(0.5 * (2 ** attempt), maximum)


async def _request_json(
    client: httpx.AsyncClient,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    max_retries: int = 5,
    max_backoff_seconds: float = 20,
) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(max_retries + 1):
        response: httpx.Response | None = None
        try:
            response = await client.get(url, params=params)
            if response.status_code not in RETRYABLE_STATUS_CODES:
                response.raise_for_status()
                payload = response.json()
                if not isinstance(payload, dict):
                    raise ChEMBLRequestError("ChEMBL returned a non-object JSON payload")
                return payload
            last_error = httpx.HTTPStatusError(
                f"retryable ChEMBL response {response.status_code}",
                request=response.request,
                response=response,
            )
        except (httpx.TransportError, httpx.TimeoutException, ValueError) as exc:
            last_error = exc

        if attempt == max_retries:
            break
        delay = _retry_delay(response, attempt, max_backoff_seconds)
        logger.warning(
            "ChEMBL request failed; retrying",
            extra={"url": url, "attempt": attempt + 1, "delay_seconds": delay},
        )
        await asyncio.sleep(delay)
    raise ChEMBLRequestError(f"ChEMBL request failed after {max_retries + 1} attempts: {url}") from last_error


async def iter_activities(
    target_chembl_id: str,
    standard_type: str = "IC50",
    *,
    page_size: int = MAX_PAGE_SIZE,
    max_records: int | None = None,
    client: httpx.AsyncClient | None = None,
    max_retries: int = 5,
) -> AsyncIterator[dict[str, Any]]:
    """Yield every activity page exposed by ChEMBL's `page_meta.next` cursor.

    `max_records` is an explicit caller safety limit. The default is unbounded and
    therefore does not silently stop at the first 1,000 records.
    """
    if page_size < 1 or page_size > MAX_PAGE_SIZE:
        raise ValueError(f"page_size must be between 1 and {MAX_PAGE_SIZE}")
    params: dict[str, Any] | None = {
        "target_chembl_id": target_chembl_id,
        "standard_type": standard_type,
        "pchembl_value__isnull": "false",
        "limit": page_size,
        "offset": 0,
        "format": "json",
    }
    url = urljoin(CHEMBL_BASE, "activity.json")
    owns_client = client is None
    active_client = client or httpx.AsyncClient(timeout=60, follow_redirects=True)
    yielded = 0
    page_number = 0
    try:
        while url:
            page_number += 1
            payload = await _request_json(active_client, url, params=params, max_retries=max_retries)
            activities = payload.get("activities")
            page_meta = payload.get("page_meta")
            if not isinstance(activities, list) or not isinstance(page_meta, dict):
                raise ChEMBLRequestError("ChEMBL activity response lacks activities/page_meta")
            logger.info("Fetched ChEMBL page", extra={"page": page_number, "records": len(activities)})
            for activity in activities:
                if not isinstance(activity, dict):
                    logger.warning("Skipping malformed ChEMBL activity", extra={"page": page_number})
                    continue
                yield activity
                yielded += 1
                if max_records is not None and yielded >= max_records:
                    return
            next_url = page_meta.get("next")
            if not next_url:
                return
            url = urljoin(CHEMBL_BASE, str(next_url))
            params = None
    finally:
        if owns_client:
            await active_client.aclose()


async def fetch_activities(
    target_chembl_id: str,
    standard_type: str = "IC50",
    limit: int | None = None,
) -> list[dict[str, Any]]:
    """Compatibility wrapper; prefer `iter_activities` for large downloads."""
    return [
        record
        async for record in iter_activities(
            target_chembl_id,
            standard_type,
            max_records=limit,
        )
    ]


async def fetch_molecule(molecule_chembl_id: str) -> dict[str, Any]:
    async with httpx.AsyncClient(timeout=30, follow_redirects=True) as client:
        return await _request_json(client, urljoin(CHEMBL_BASE, f"molecule/{molecule_chembl_id}.json"))


async def iter_assays_by_ids(
    assay_chembl_ids: Sequence[str],
    *,
    batch_size: int = 50,
    client: httpx.AsyncClient | None = None,
    max_retries: int = 5,
) -> AsyncIterator[tuple[dict[str, Any], str]]:
    """Yield ChEMBL assay records in bounded batches with their request URL.

    The returned URL is retained by the caller as record-level provenance. IDs are
    sorted and deduplicated so repeated runs produce the same request sequence.
    """
    if batch_size < 1 or batch_size > 100:
        raise ValueError("batch_size must be between 1 and 100")
    identifiers = sorted({str(value).strip() for value in assay_chembl_ids if str(value).strip()})
    owns_client = client is None
    active_client = client or httpx.AsyncClient(timeout=60, follow_redirects=True)
    base_url = urljoin(CHEMBL_BASE, "assay.json")
    try:
        for start in range(0, len(identifiers), batch_size):
            batch = identifiers[start:start + batch_size]
            params: dict[str, Any] | None = {
                "assay_chembl_id__in": ",".join(batch),
                "limit": batch_size,
                "offset": 0,
                "format": "json",
            }
            url = base_url
            while url:
                payload = await _request_json(active_client, url, params=params, max_retries=max_retries)
                assays = payload.get("assays")
                page_meta = payload.get("page_meta")
                if not isinstance(assays, list) or not isinstance(page_meta, dict):
                    raise ChEMBLRequestError("ChEMBL assay response lacks assays/page_meta")
                request_url = str(active_client.build_request("GET", url, params=params).url)
                for assay in assays:
                    if isinstance(assay, dict):
                        yield assay, request_url
                    else:
                        logger.warning("Skipping malformed ChEMBL assay record")
                next_url = page_meta.get("next")
                url = urljoin(CHEMBL_BASE, str(next_url)) if next_url else ""
                params = None
    finally:
        if owns_client:
            await active_client.aclose()

