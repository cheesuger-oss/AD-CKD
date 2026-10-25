#!/usr/bin/env python3
"""Freeze six GO-based immune activation programs from the local OrgDb."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sqlite3
from pathlib import Path


PROGRAMS = [
    ("IMM-MAC", "macrophage activation", "GO:0042116"),
    ("IMM-NEU", "neutrophil activation", "GO:0042119"),
    ("IMM-T", "T cell activation", "GO:0042110"),
    ("IMM-B", "B cell activation", "GO:0042113"),
    ("IMM-NK", "natural killer cell activation", "GO:0030101"),
    ("IMM-DC", "dendritic cell activation", "GO:0036336"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def freeze_programs(db_path: Path, output_path: Path, audit_path: Path):
    db_path = Path(db_path)
    db_hash = sha256(db_path)
    connection = sqlite3.connect(db_path)
    metadata = dict(connection.execute("select name, value from metadata"))
    rows: list[dict[str, str]] = []
    summary: list[dict[str, object]] = []
    for program_id, program_name, go_id in PROGRAMS:
        genes = connection.execute(
            "select distinct i.symbol, g.gene_id "
            "from go_bp_all x join gene_info i using(_id) join genes g using(_id) "
            "where x.go_id=? and i.symbol is not null and i.symbol<>'' order by i.symbol",
            (go_id,),
        ).fetchall()
        if len(genes) < 10:
            raise ValueError(f"immune program too small or absent: {go_id}")
        summary.append({"program_id": program_id, "program_name": program_name, "go_id": go_id, "gene_n": len(genes)})
        for symbol, entrez in genes:
            rows.append({
                "program_id": program_id,
                "program_name": program_name,
                "go_id": go_id,
                "human_gene": symbol,
                "entrez_id": entrez,
                "ontology_relation": "GOALL (term plus descendants)",
                "interpretation": "activation_program_not_cell_fraction",
                "source_database": "org.Hs.eg.db 3.22.0 local SQLite",
                "go_source_date": metadata.get("GOSOURCEDATE", "UNKNOWN"),
                "entrez_go_source_date": metadata.get("GOEGSOURCEDATE", "UNKNOWN"),
                "source_sha256": db_hash,
            })
    connection.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    audit_path.write_text(
        "# Phase 4 immune-program freeze audit\n\n"
        "Frozen: 2026-08-09\n\n"
        f"- Source: local `org.Hs.eg.db` 3.22.0 SQLite; SHA-256 `{db_hash}`.\n"
        f"- Gene Ontology source date: {metadata.get('GOSOURCEDATE', 'UNKNOWN')}.\n"
        f"- Entrez-to-GO source date: {metadata.get('GOEGSOURCEDATE', 'UNKNOWN')}.\n"
        "- Membership uses `go_bp_all`, therefore each program includes the named term and descendants.\n"
        "- Planned sample score: mean within-sample percentile rank of the observed program genes (singscore-style rank score), computed independently within each dataset.\n"
        "- Interpretation ceiling: transcriptional activation program; not a direct immune-cell fraction or histological cell count.\n\n"
        "| Program | GO ID | Frozen genes |\n|---|---|---:|\n" +
        "\n".join(f"| {x['program_name']} | {x['go_id']} | {x['gene_n']} |" for x in summary) +
        "\n\nAll six programs form one within-dataset BH correction family and must be reported regardless of direction or significance.\n",
        encoding="utf-8",
    )
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = freeze_programs(args.db, args.out, args.audit)
    print(f"frozen immune programs: {len(summary)}; membership rows: {len(rows)}")


if __name__ == "__main__":
    main()
