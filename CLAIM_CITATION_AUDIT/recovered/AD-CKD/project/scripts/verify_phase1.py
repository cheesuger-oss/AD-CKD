#!/usr/bin/env python3
"""Deterministic Phase 1 artifact and result gate."""

from __future__ import annotations

import csv
import hashlib
import sys
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "project" / "results" / "geo"


def read_csv(name: str, delimiter: str = ",") -> list[dict[str, str]]:
    path = OUT / name
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    checks: list[tuple[str, bool, str]] = []

    required = [
        "GSE146638_DEG.csv", "GSE325608_DEG.csv", "intersection_genes.csv",
        "convergence_pathways.csv", "ortholog_rat_human_1to1.csv",
        "ortholog_mouse_human_1to1.csv", "PHASE1_METADATA_AUDIT.md", "PHASE1_RESULTS.md",
        "figures/GSE146638_PCA.png", "figures/GSE325608_PCA.png",
        "figures/GSE146638_volcano.png", "figures/GSE325608_volcano.png",
        "figures/convergence_pathway_heatmap.png",
    ]
    for name in required:
        path = OUT / name
        checks.append((f"artifact:{name}", path.exists() and path.stat().st_size > 100, str(path.stat().st_size if path.exists() else 0)))

    m1 = read_csv("GSE146638_sample_metadata.tsv", "\t")
    m2 = read_csv("GSE325608_sample_metadata.tsv", "\t")
    checks.extend([
        ("GSE146638_samples_14", len(m1) == 14, str(len(m1))),
        ("GSE146638_primary_5v5", sum(x["primary_group"] == "Control" and x["include_primary"] == "yes" for x in m1) == 5 and sum(x["primary_group"] == "Uremic" and x["include_primary"] == "yes" for x in m1) == 5, "frozen"),
        ("GSE325608_samples_48", len(m2) == 48, str(len(m2))),
        ("GSE325608_primary_20v16", sum(x["primary_group"] == "Control" and x["include_primary"] == "yes" for x in m2) == 20 and sum(x["primary_group"] == "AAD" and x["include_primary"] == "yes" for x in m2) == 16, "frozen"),
        ("GSE325608_nonAD_12", sum(x["primary_group"] == "NonAD" for x in m2) == 12, str(sum(x["primary_group"] == "NonAD" for x in m2))),
        ("GSE325608_all_GSM_unique", len({x["geo_accession"] for x in m2}) == 48, str(len({x["geo_accession"] for x in m2}))),
    ])

    rat = read_csv("ortholog_rat_human_1to1.csv")
    mouse = read_csv("ortholog_mouse_human_1to1.csv")
    for label, mapping in (("rat", rat), ("mouse", mouse)):
        valid = all(x["orthology_type"] == "ortholog_one2one" and float(x["orthology_confidence"]) == 1 for x in mapping)
        checks.append((f"{label}_mapping_strict_1to1", valid and len(mapping) > 1000, str(len(mapping))))

    metrics = {x["metric"]: x["value"] for x in read_csv("phase1_convergence_summary.tsv", "\t")}
    shared = int(metrics["shared_FDR_lt_0.01_pathways"])
    focus = int(metrics["focus_pathways"])
    checks.extend([
        ("shared_pathways_at_least_5", shared >= 5, str(shared)),
        ("focus_pathways_at_least_3", focus >= 3, str(focus)),
        ("support_gate", metrics["support"] == "SUPPORTED", metrics["support"]),
    ])
    kegg_verified = metrics["KEGG_status"].startswith("VERIFIED")
    checks.append(("KEGG_annotation_verified", kegg_verified, metrics["KEGG_status"]))

    failures = [x for x in checks if not x[1]]
    core_failures = [x for x in failures if x[0] != "KEGG_annotation_verified"]
    status = "PASS" if not failures else ("PASS_WITH_KEGG_PENDING" if not core_failures else "FAIL")
    hashes = []
    for name in required:
        path = OUT / name
        if path.exists():
            hashes.append((name, sha256(path)))

    lines = [
        "# Phase 1 验证状态", "", f"日期：{date.today().isoformat()}", "", f"状态：**{status}**", "",
        "## 检查", "",
    ]
    lines += [f"- [{'x' if ok else ' '}] {name}: {detail}" for name, ok, detail in checks]
    lines += ["", "## 核心产物 SHA-256", ""]
    lines += [f"- `{name}`: `{digest}`" for name, digest in hashes]
    (OUT / "PHASE1_STATUS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Phase 1 verification: {status}; checks={len(checks)}; failures={len(failures)}")
    for name, _, detail in failures:
        print(f"FAIL\t{name}\t{detail}")
    return 0 if status == "PASS" else (2 if status == "PASS_WITH_KEGG_PENDING" else 1)


if __name__ == "__main__":
    sys.exit(main())
