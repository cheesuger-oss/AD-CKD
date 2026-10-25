#!/usr/bin/env python3
"""Freeze clumped renal instruments and prepare one union outcome extraction list."""

from __future__ import annotations

import csv
import json
from pathlib import Path


EXPOSURES = (
    "serum_cystatin_c", "egfrcys", "bun_ea", "ckd_ea", "microalbuminuria",
    "rapid3", "ckdi25", "uacr_ea",
)


def clumped_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if not lines:
        return []
    header = lines[0]
    snp_idx = header.index("SNP")
    return [fields[snp_idx] for fields in lines[1:]]


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    base = root / "project/data/gwas/processed/renal_family"
    union: set[str] = set()
    summary = []
    for exposure in EXPOSURES:
        ids = clumped_ids(base / "clump" / f"{exposure}.clumped")
        union.update(ids)
        with (base / f"{exposure}_gws_candidates.tsv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        by_snp = {row["SNP"]: row for row in rows}
        selected = [by_snp[snp] for snp in ids]
        with (base / f"{exposure}_clumped_instruments.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected)
        summary.append({"exposure_id": exposure, "clumped_instruments": len(ids)})
    (base / "renal_family_union_rsids.txt").write_text("".join(f"{x}\n" for x in sorted(union)), encoding="utf-8")
    payload = {"exposures": summary, "unique_union_rsids": len(union)}
    (base / "clump_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
