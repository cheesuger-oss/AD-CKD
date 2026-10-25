#!/usr/bin/env python3
"""Audit locked local external AAD cohorts before any model fitting."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


COHORTS = [
    ("GSE52093", "GSE52093_series_matrix.txt.gz", "GPL10558", 1),
    ("GSE98770", "GSE98770-GPL14550_series_matrix.txt.gz", "GPL14550", 2),
    ("GSE147026", "GSE147026_series_matrix.txt.gz", "GPL24676", 3),
    ("GSE190635", "GSE190635_series_matrix.txt.gz", "GPL570", 4),
]

FIELDS = [
    "execution_order", "dataset_id", "platform", "assay", "case_n", "control_n",
    "sample_n", "expression_rows", "matrix_sha256", "family_archive_sha256",
    "metadata_label_source", "annotation_source", "fitness", "limitation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_matrix(path: Path) -> dict[str, object]:
    metadata: dict[str, list[str]] = {}
    sample_ids: list[str] = []
    expression_rows = 0
    inside = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "!series_matrix_table_begin":
                inside = True
                continue
            if line == "!series_matrix_table_end":
                inside = False
                continue
            if inside:
                if expression_rows == 0:
                    sample_ids = [value.strip('"') for value in line.split("\t")[1:]]
                expression_rows += 1
                continue
            if line.startswith("!Sample_"):
                parts = line.split("\t")
                metadata.setdefault(parts[0], []).append([value.strip('"') for value in parts[1:]])
    if expression_rows:
        expression_rows -= 1  # table header
    accessions = metadata.get("!Sample_geo_accession", [[]])[0]
    if sample_ids and sample_ids != accessions:
        raise ValueError(f"matrix/sample metadata order mismatch: {path.name}")
    return {"metadata": metadata, "sample_ids": accessions, "expression_rows": expression_rows}


def classify_groups(dataset_id: str, metadata: dict[str, list[list[str]]]) -> tuple[int, int]:
    titles = metadata.get("!Sample_title", [[]])[0]
    characteristics = metadata.get("!Sample_characteristics_ch1", [])
    labels: list[str] = []
    for index, title in enumerate(titles):
        text = " ".join([title] + [row[index] for row in characteristics if index < len(row)]).lower()
        if any(token in text for token in ["aortic dissection", "aorta dissected", "ataad"]):
            labels.append("case")
        elif any(token in text for token in ["normal", "healthy", "donor", "control"]):
            labels.append("control")
        else:
            raise ValueError(f"unclassified sample in {dataset_id}: {title}")
    return labels.count("case"), labels.count("control")


def audit_external_cohorts(raw_dir: Path, output_tsv: Path, output_md: Path) -> list[dict[str, str]]:
    raw_dir = Path(raw_dir)
    rows: list[dict[str, str]] = []
    limitations = {
        "GSE52093": "small surgical case-control cohort; one female case; control selection differs",
        "GSE98770": "small cohort; transplant-donor controls; mRNA platform only",
        "GSE147026": "series matrix contains metadata only and zero expression rows; processed expression files absent locally",
        "GSE190635": "very small all-male cohort; Affymetrix platform and control selection differ",
    }
    assays = {
        "GSE52093": "Illumina HumanHT-12 v4 expression array",
        "GSE98770": "Agilent SurePrint G3 human GE 8x60K mRNA array",
        "GSE147026": "bulk RNA-seq metadata",
        "GSE190635": "Affymetrix HG-U133 Plus 2.0 array",
    }
    for dataset_id, filename, platform, order in COHORTS:
        matrix = raw_dir / filename
        family = raw_dir / f"{dataset_id}_family.xml.tgz"
        if not matrix.exists() or not family.exists():
            raise FileNotFoundError(f"missing locked input for {dataset_id}")
        parsed = parse_matrix(matrix)
        case_n, control_n = classify_groups(dataset_id, parsed["metadata"])
        sample_n = len(parsed["sample_ids"])
        if case_n + control_n != sample_n or len(set(parsed["sample_ids"])) != sample_n:
            raise ValueError(f"label or duplicate-sample failure: {dataset_id}")
        expression_rows = int(parsed["expression_rows"])
        rows.append({
            "execution_order": str(order),
            "dataset_id": dataset_id,
            "platform": platform,
            "assay": assays[dataset_id],
            "case_n": str(case_n),
            "control_n": str(control_n),
            "sample_n": str(sample_n),
            "expression_rows": str(expression_rows),
            "matrix_sha256": sha256(matrix),
            "family_archive_sha256": sha256(family),
            "metadata_label_source": "local GEO series matrix Sample_title/characteristics",
            "annotation_source": f"local {dataset_id}_family.xml.tgz platform table",
            "fitness": "READY_WITH_LIMITATIONS" if expression_rows > 0 else "NOT_READY",
            "limitation": limitations[dataset_id],
        })

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    output_md.write_text(
        "# Phase 4 external AAD cohort audit\n\n"
        "Frozen: 2026-08-09\n\n"
        "The execution order was fixed before reading model outcomes. Local series matrices and family archives were hashed, and labels were derived from GEO sample titles/characteristics.\n\n"
        "| Order | Dataset | Platform | Cases | Controls | Expression rows | Fitness |\n"
        "|---:|---|---|---:|---:|---:|---|\n" +
        "\n".join(
            f"| {r['execution_order']} | {r['dataset_id']} | {r['platform']} | {r['case_n']} | {r['control_n']} | {r['expression_rows']} | {r['fitness']} |"
            for r in rows
        ) +
        "\n\n## Decisions\n\n"
        "- GSE52093, GSE98770 and GSE190635 are locally executable with strong small-sample and control-selection limitations.\n"
        "- GSE147026 is retained in the complete audit but is not executable from the current metadata-only matrix; it cannot be silently dropped or called a negative validation.\n"
        "- No cohort has been evaluated for model performance at this audit stage.\n"
        "- External validation cannot upgrade a retrospective tissue classifier to incident-risk prediction.\n",
        encoding="utf-8",
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows = audit_external_cohorts(args.raw_dir, args.out, args.report)
    print(f"audited external AAD cohorts: {len(rows)}")


if __name__ == "__main__":
    main()
