#!/usr/bin/env python3
"""Freeze the renal-only Phase 4 candidate universe without AAD labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


FIELDS = [
    "candidate_id",
    "candidate_origin",
    "aad_label_used",
    "renal_dataset",
    "rat_gene",
    "human_gene",
    "renal_log_fc",
    "renal_direction",
    "renal_adj_p",
    "orthology_type",
    "orthology_confidence",
    "source_symbol_mapping_count",
    "mapping_note",
    "selection_rule",
    "deg_input_sha256",
    "ortholog_input_sha256",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def freeze_candidates(deg_path: Path, ortholog_path: Path, output_path: Path, audit_path: Path) -> list[dict[str, str]]:
    deg_path = Path(deg_path)
    ortholog_path = Path(ortholog_path)
    output_path = Path(output_path)
    audit_path = Path(audit_path)

    deg_rows = _read_csv(deg_path)
    map_rows = _read_csv(ortholog_path)
    required_deg = {"gene_symbol", "logFC", "adj.P.Val"}
    required_map = {"source_symbol", "human_symbol", "orthology_type", "orthology_confidence"}
    if not deg_rows or not required_deg.issubset(deg_rows[0]):
        raise ValueError(f"DEG input lacks required fields: {sorted(required_deg)}")
    if not map_rows or not required_map.issubset(map_rows[0]):
        raise ValueError(f"ortholog input lacks required fields: {sorted(required_map)}")

    selected_deg: dict[str, dict[str, str]] = {}
    for row in deg_rows:
        gene = row["gene_symbol"].strip()
        log_fc = float(row["logFC"])
        adj_p = float(row["adj.P.Val"])
        if adj_p < 0.05 and abs(log_fc) > 0.5:
            if gene in selected_deg:
                raise ValueError(f"duplicate DEG gene symbol: {gene}")
            selected_deg[gene] = row

    valid_maps = [
        row for row in map_rows
        if row["source_symbol"].strip() in selected_deg
        and row["human_symbol"].strip()
        and row["orthology_type"].strip() == "ortholog_one2one"
        and row["orthology_confidence"].strip() == "1"
    ]
    source_counts = Counter(row["source_symbol"].strip() for row in valid_maps)
    human_sources: dict[str, set[str]] = {}
    for row in valid_maps:
        human_sources.setdefault(row["human_symbol"].strip(), set()).add(row["source_symbol"].strip())
    conflicts = {human: sources for human, sources in human_sources.items() if len(sources) > 1}
    if conflicts:
        first = sorted(conflicts.items())[0]
        raise ValueError(f"human gene maps from multiple renal source genes: {first[0]} <- {sorted(first[1])}")

    deg_hash = sha256(deg_path)
    map_hash = sha256(ortholog_path)
    rows: list[dict[str, str]] = []
    seen_human: set[str] = set()
    for mapping in sorted(valid_maps, key=lambda r: (r["human_symbol"].strip(), r["source_symbol"].strip())):
        rat_gene = mapping["source_symbol"].strip()
        human_gene = mapping["human_symbol"].strip()
        if human_gene in seen_human:
            raise ValueError(f"duplicate human candidate after sorting: {human_gene}")
        seen_human.add(human_gene)
        deg = selected_deg[rat_gene]
        log_fc = float(deg["logFC"])
        count = source_counts[rat_gene]
        rows.append({
            "candidate_id": f"RENAL-{len(rows) + 1:04d}",
            "candidate_origin": "renal_only",
            "aad_label_used": "false",
            "renal_dataset": "GSE146638",
            "rat_gene": rat_gene,
            "human_gene": human_gene,
            "renal_log_fc": f"{log_fc:.12g}",
            "renal_direction": "up" if log_fc > 0 else "down",
            "renal_adj_p": f"{float(deg['adj.P.Val']):.12g}",
            "orthology_type": mapping["orthology_type"].strip(),
            "orthology_confidence": mapping["orthology_confidence"].strip(),
            "source_symbol_mapping_count": str(count),
            "mapping_note": "source symbol appears in multiple one-to-one rows" if count > 1 else "unique source symbol row",
            "selection_rule": "GSE146638 adj.P.Val<0.05 & abs(logFC)>0.5; human ortholog_one2one confidence=1",
            "deg_input_sha256": deg_hash,
            "ortholog_input_sha256": map_hash,
        })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    multi_source_symbols = sum(1 for count in source_counts.values() if count > 1)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(
        "\n".join([
            "# Phase 4 renal-origin candidate freeze audit",
            "",
            "Frozen: 2026-08-09",
            "",
            f"- DEG input: `{deg_path}`; SHA-256 `{deg_hash}`.",
            f"- Orthology input: `{ortholog_path}`; SHA-256 `{map_hash}`.",
            f"- Renal DEGs passing the frozen rule: {len(selected_deg)}.",
            f"- Valid one-to-one mapping rows: {len(valid_maps)}.",
            f"- Unique frozen human candidates: {len(rows)}.",
            f"- Rat source symbols represented by multiple one-to-one mapping rows: {multi_source_symbols}; retained with an explicit row-level flag because each human target remains unique.",
            "- AAD labels or AAD expression values used in candidate construction: **no**.",
            "- Candidate role: primary renal-origin feature universe for downstream fold-internal modeling; not a causal gene set.",
            "",
            "## Fail-closed checks",
            "",
            "- Duplicate renal DEG symbols: rejected.",
            "- One human symbol reached from multiple renal source genes: rejected.",
            "- Non-one-to-one or confidence<1 mapping rows: excluded.",
            "- Threshold drift from adj.P.Val<0.05 and |logFC|>0.5: rejected by tests.",
            "",
        ]) + "\n",
        encoding="utf-8",
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deg", type=Path, required=True)
    parser.add_argument("--ortholog", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    rows = freeze_candidates(args.deg, args.ortholog, args.out, args.audit)
    print(f"frozen renal-only candidates: {len(rows)}")


if __name__ == "__main__":
    main()
