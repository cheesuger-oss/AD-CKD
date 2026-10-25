#!/usr/bin/env python3
"""Query DGIdb for the frozen top-20 descriptive topology nodes."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import ssl
import urllib.request
from pathlib import Path


ENDPOINT = "https://dgidb.org/api/graphql"
TLS_CONTEXT = ssl.create_default_context(cafile="/etc/ssl/cert.pem")
TOP_N = 20


def flatten(payload: dict, queried: list[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    nodes = payload.get("data", {}).get("genes", {}).get("nodes", [])
    by_name = {str(node.get("name", "")).upper(): node for node in nodes}
    audit, interactions = [], []
    for rank, gene in enumerate(queried, 1):
        node = by_name.get(gene.upper())
        items = [] if node is None else node.get("interactions", []) or []
        audit.append({"topology_rank": rank, "query_gene": gene,
                      "dgidb_mapping_status": "mapped" if node else "unmapped",
                      "interaction_count": len(items)})
        for item in items:
            types = item.get("interactionTypes", []) or []
            interactions.append({
                "topology_rank": rank, "gene": gene,
                "drug": (item.get("drug") or {}).get("name", ""),
                "drug_concept_id": (item.get("drug") or {}).get("conceptId", ""),
                "interaction_score": item.get("interactionScore", ""),
                "interaction_types": "|".join(sorted({str(x.get("type", "")) for x in types if x.get("type")})),
                "directionality": "|".join(sorted({str(x.get("directionality", "")) for x in types if x.get("directionality")})),
                "sources": "|".join(sorted({str(x.get("sourceDbName", "")) for x in item.get("sources", []) if x.get("sourceDbName")})),
                "pmids": "|".join(sorted({str(x.get("pmid", "")) for x in item.get("publications", []) if x.get("pmid")})),
            })
    interactions.sort(key=lambda row: (int(row["topology_rank"]), -float(row["interaction_score"] or 0), str(row["drug"])))
    return audit, interactions


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    base = root / "project/results/network"
    with (base / "hub_genes.csv").open(newline="", encoding="utf-8") as handle:
        queried = [row["gene"] for row in list(csv.DictReader(handle))[:TOP_N]]
    names = ",".join(json.dumps(name) for name in queried)
    query = """{ genes(names: [%s]) { nodes { name conceptId interactions { drug { name conceptId } interactionScore interactionTypes { type directionality } publications { pmid } sources { sourceDbName } } } } }""" % names
    body = json.dumps({"query": query}).encode("utf-8")
    request = urllib.request.Request(ENDPOINT, data=body, method="POST",
        headers={"Content-Type": "application/json", "User-Agent": "AD-CKD/phase3"})
    with urllib.request.urlopen(request, timeout=240, context=TLS_CONTEXT) as response:
        raw = response.read()
    payload = json.loads(raw)
    if payload.get("errors"):
        raise ValueError(payload["errors"])
    audit, interactions = flatten(payload, queried)
    (base / "dgidb_raw.json").write_bytes(raw)
    for name, rows in (("dgidb_gene_audit.tsv", audit), ("drug_targets.csv", interactions)):
        path = base / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(),
                                    delimiter="\t" if name.endswith(".tsv") else ",", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    metadata = {
        "schema": "ad-ckd.phase3.dgidb.v1", "queried_at": dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds"),
        "endpoint": ENDPOINT, "reported_site_version": "5.0.11 (2026-01-13)",
        "selection_rule": "first 20 nodes from frozen degree-descending topology table",
        "queried_genes": queried, "mapped_gene_count": sum(x["dgidb_mapping_status"] == "mapped" for x in audit),
        "interaction_row_count": len(interactions), "raw_sha256": hashlib.sha256(raw).hexdigest(),
        "interpretation": "database interactions are hypothesis-generating; score is not treatment efficacy or clinical recommendation",
    }
    (base / "dgidb_query_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
