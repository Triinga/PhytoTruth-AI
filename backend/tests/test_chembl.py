import asyncio

import httpx

from app import chembl


def test_iter_activities_follows_every_next_page():
    requests = []

    def handler(request: httpx.Request):
        requests.append(str(request.url))
        if "offset=2" in str(request.url):
            return httpx.Response(200, json={"activities": [{"activity_id": 3}], "page_meta": {"next": None}})
        return httpx.Response(
            200,
            json={
                "activities": [{"activity_id": 1}, {"activity_id": 2}],
                "page_meta": {"next": "/chembl/api/data/activity.json?limit=2&offset=2"},
            },
        )

    async def collect():
        async with httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False) as client:
            return [item async for item in chembl.iter_activities("CHEMBL335", page_size=2, client=client)]

    records = asyncio.run(collect())
    assert [record["activity_id"] for record in records] == [1, 2, 3]
    assert len(requests) == 2
    assert "offset=2" in requests[1]


def test_iter_activities_retries_429(monkeypatch):
    attempts = 0
    sleeps = []

    def handler(request: httpx.Request):
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            return httpx.Response(429, headers={"Retry-After": "0"})
        return httpx.Response(200, json={"activities": [{"activity_id": 7}], "page_meta": {"next": None}})

    async def fake_sleep(delay):
        sleeps.append(delay)

    monkeypatch.setattr(chembl.asyncio, "sleep", fake_sleep)

    async def collect():
        async with httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False) as client:
            return [item async for item in chembl.iter_activities("CHEMBL335", client=client)]

    records = asyncio.run(collect())
    assert records[0]["activity_id"] == 7
    assert attempts == 2
    assert sleeps == [0.0]


def test_explicit_max_records_is_honored():
    def handler(request: httpx.Request):
        return httpx.Response(
            200,
            json={"activities": [{"activity_id": 1}, {"activity_id": 2}], "page_meta": {"next": "/next"}},
        )

    async def collect():
        async with httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False) as client:
            return [item async for item in chembl.iter_activities("CHEMBL335", max_records=1, client=client)]

    assert len(asyncio.run(collect())) == 1


def test_iter_assays_batches_deduplicates_and_preserves_request_url():
    requests = []

    def handler(request: httpx.Request):
        requests.append(str(request.url))
        ids = request.url.params["assay_chembl_id__in"].split(",")
        return httpx.Response(
            200,
            json={"assays": [{"assay_chembl_id": value} for value in ids], "page_meta": {"next": None}},
        )

    async def collect():
        async with httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False) as client:
            return [
                item async for item in chembl.iter_assays_by_ids(
                    ["CHEMBL3", "CHEMBL1", "CHEMBL1", "CHEMBL2"], batch_size=2, client=client
                )
            ]

    records = asyncio.run(collect())
    assert [record[0]["assay_chembl_id"] for record in records] == ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    assert len(requests) == 2
    assert all(record_url.startswith("https://www.ebi.ac.uk/chembl/api/data/assay.json?") for _, record_url in records)

