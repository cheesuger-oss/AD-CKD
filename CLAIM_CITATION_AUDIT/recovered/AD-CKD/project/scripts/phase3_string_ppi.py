#!/usr/bin/env python3
"""Map frozen genes and retrieve a version-pinned STRING functional network."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import ssl
import urllib.parse
import urllib.request
from pathlib import Path


SPECIES = 9606
REQUIRED_SCORE = 400
NETWORK_TYPE = "functional"
CALLER_IDENTITY = "AD-CKD-phase3-local-research"
TLS_CONTEXT = ssl.create_default_context(cafile="/etc/ssl/cert.pem")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def post(url: str, params: dict[str, object], timeout: int = 240) -> bytes:
    body = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=body,
        headers={"User-Agent": "AD-CKD/phase3", "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout, context=TLS_CONTEXT) as response:
        return response.read()


def get(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "AD-CKD/phase3"})
    with urllib.request.urlopen(request, timeout=timeout, context=TLS_CONTEXT) as response:
        return response.read()


def parse_tsv(data: bytes) -> tuple[list[str], list[dict[str, str]]]:
    text = data.decode("utf-8-sig")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("STRING returned an empty or headerless TSV response")
    return reader.fieldnames, list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(root: Path) -> None:
    output = root / "project/results/network"
    output.mkdir(parents=True, exist_ok=True)
    gene_path = output / "union_genes.txt"
    genes = gene_path.read_text(encoding="utf-8").splitlines()
    if not genes or genes != sorted(set(genes)):
        raise ValueError("union_genes.txt must be non-empty, unique, and sorted")

    queried_at = dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")
    version_raw = get("https://string-db.org/api/json/version")
    version_json = json.loads(version_raw)
    if not isinstance(version_json, list) or len(version_json) != 1:
        raise ValueError(f"unexpected STRING version response: {version_json!r}")
    version = version_json[0]
    stable_value = version.get("stable_address") or version.get("string_stable_address")
    if not stable_value:
        raise ValueError(f"STRING version response lacks a stable address: {version!r}")
    stable = str(stable_value).rstrip("/")
    api = stable + "/api"

    mapping_params = {
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    mapping_raw = post(api + "/tsv/get_string_ids", mapping_params)
    mapping_fields, mapping_rows = parse_tsv(mapping_raw)
    required_mapping = {"queryItem", "queryIndex", "stringId", "ncbiTaxonId", "preferredName"}
    if not required_mapping.issubset(mapping_fields):
        raise ValueError(f"mapping fields changed: {mapping_fields}")

    by_query: dict[str, list[dict[str, str]]] = {gene: [] for gene in genes}
    for row in mapping_rows:
        query = row["queryItem"]
        if query in by_query:
            by_query[query].append(row)
    audit_rows: list[dict[str, object]] = []
    mapped_ids: list[str] = []
    for gene in genes:
        matches = by_query[gene]
        mapped_ids.extend(row["stringId"] for row in matches)
        audit_rows.append({
            "query_gene": gene,
            "mapping_status": "mapped" if len(matches) == 1 else ("unmapped" if not matches else "multiple"),
            "mapping_count": len(matches),
            "string_id": "|".join(row["stringId"] for row in matches),
            "preferred_name": "|".join(row["preferredName"] for row in matches),
            "ncbi_taxon_id": "|".join(row["ncbiTaxonId"] for row in matches),
        })
    mapped_ids = sorted(set(mapped_ids))
    if not mapped_ids:
        raise ValueError("STRING mapped zero input genes")

    network_params = {
        "identifiers": "\r".join(mapped_ids),
        "species": SPECIES,
        "required_score": REQUIRED_SCORE,
        "network_type": NETWORK_TYPE,
        "caller_identity": CALLER_IDENTITY,
    }
    network_raw = post(api + "/tsv/network", network_params)
    network_fields, network_rows = parse_tsv(network_raw)
    required_network = {
        "stringId_A", "stringId_B", "preferredName_A", "preferredName_B",
        "ncbiTaxonId", "score",
    }
    if not required_network.issubset(network_fields):
        raise ValueError(f"network fields changed: {network_fields}")

    edges: dict[tuple[str, str], float] = {}
    for row in network_rows:
        if int(row["ncbiTaxonId"]) != SPECIES:
            raise ValueError("network response contains a non-human edge")
        score = float(row["score"])
        if score + 1e-12 < REQUIRED_SCORE / 1000:
            raise ValueError(f"network response contains score below threshold: {score}")
        a, b = sorted((row["preferredName_A"], row["preferredName_B"]))
        if a == b:
            continue
        edges[(a, b)] = max(score, edges.get((a, b), -1.0))

    (output / "string_version_raw.json").write_bytes(version_raw)
    (output / "string_mapping_raw.tsv").write_bytes(mapping_raw)
    (output / "string_network_raw.tsv").write_bytes(network_raw)
    write_tsv(
        output / "string_id_mapping.tsv",
        ["query_gene", "mapping_status", "mapping_count", "string_id", "preferred_name", "ncbi_taxon_id"],
        audit_rows,
    )
    write_tsv(
        output / "string_unmapped_genes.tsv",
        ["query_gene", "mapping_status", "mapping_count", "string_id", "preferred_name", "ncbi_taxon_id"],
        [row for row in audit_rows if row["mapping_status"] != "mapped"],
    )
    with (output / "union_ppi.sif").open("w", encoding="utf-8", newline="") as handle:
        handle.write("source\ttarget\tscore\n")
        for (a, b), score in sorted(edges.items()):
            handle.write(f"{a}\t{b}\t{score:.3f}\n")

    mapped_names = {row["preferredName"] for row in mapping_rows}
    connected_names = {name for edge in edges for name in edge}
    metadata = {
        "schema": "ad-ckd.phase3.string.v1",
        "queried_at": queried_at,
        "string_version": version["string_version"],
        "string_stable_address": stable,
        "species": SPECIES,
        "network_type": NETWORK_TYPE,
        "required_score_api_units": REQUIRED_SCORE,
        "required_score_output_scale": REQUIRED_SCORE / 1000,
        "caller_identity": CALLER_IDENTITY,
        "input_gene_count": len(genes),
        "mapped_query_count": sum(row["mapping_status"] == "mapped" for row in audit_rows),
        "unmapped_query_count": sum(row["mapping_status"] == "unmapped" for row in audit_rows),
        "multiple_mapping_query_count": sum(row["mapping_status"] == "multiple" for row in audit_rows),
        "unique_string_id_count": len(mapped_ids),
        "mapped_preferred_name_count": len(mapped_names),
        "connected_node_count": len(connected_names),
        "isolated_mapped_node_count": len(mapped_names - connected_names),
        "edge_count": len(edges),
        "interpretation": "functional association network; topology is descriptive and not causal evidence",
        "input_sha256": sha256_bytes(gene_path.read_bytes()),
        "raw_response_sha256": {
            "version": sha256_bytes(version_raw),
            "mapping": sha256_bytes(mapping_raw),
            "network": sha256_bytes(network_raw),
        },
    }
    (output / "string_query_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    run(args.root.resolve())


if __name__ == "__main__":
    main()
