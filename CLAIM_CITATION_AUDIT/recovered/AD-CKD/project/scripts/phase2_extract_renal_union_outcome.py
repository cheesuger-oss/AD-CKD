#!/usr/bin/env python3
"""Extract the frozen renal-family rsID union from the local FinnGen outcome."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    base = root / "project/data/gwas/processed/renal_family"
    wanted = set((base / "renal_family_union_rsids.txt").read_text(encoding="utf-8").split())
    source = Path("inputs/finngen_R13_I9_AORTDIS.gz")
    output = base / "finngen_aortdis_union_rows.tsv"
    scanned = 0
    matched: dict[str, dict[str, str]] = {}
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "rsids" not in reader.fieldnames:
            raise ValueError(f"FinnGen rsids field missing: {reader.fieldnames}")
        for row in reader:
            scanned += 1
            rsids = row["rsids"].split(",") if row["rsids"] else []
            hits = wanted.intersection(rsids)
            for rsid in hits:
                record = dict(row)
                record["rsids"] = rsid
                matched.setdefault(rsid, record)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(matched[key] for key in sorted(matched))
    payload = {
        "source": str(source), "rows_scanned": scanned, "requested_unique_rsids": len(wanted),
        "matched_unique_rsids": len(matched), "missing_unique_rsids": len(wanted - matched.keys()),
        "output": str(output),
    }
    (base / "outcome_extract_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
