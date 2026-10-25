#!/usr/bin/env python3
"""Convert PLINK clumps into the frozen pairwise colocalization window family."""

import csv
import json
from pathlib import Path

from prepare_coloc_loci import merge_windows


RENAL = ["egfrcrea", "egfrcys", "cystatin_c", "bun", "ckd", "uacr"]


def read_clumps(path, trait):
    rows = []
    with open(path) as handle:
        for line in handle:
            values = line.split()
            if not values or values[0] == "CHR":
                continue
            rows.append({"trait": trait, "chrom": values[0], "pos_b37": int(values[3]), "rsid": values[2], "p": values[4]})
    return rows


def main():
    root = Path("project/results/shared_genetics/coloc")
    leads = {trait: read_clumps(root / f"{trait}_clump.clumped", trait) for trait in ["aad", *RENAL]}
    with open(root / "independent_leads.tsv", "w", newline="") as handle:
        fields = ["trait", "chrom", "pos_b37", "rsid", "p"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for trait in leads:
            writer.writerows(leads[trait])

    memberships = []
    windows = []
    for renal in RENAL:
        union = [(x["chrom"], x["pos_b37"], f'{x["trait"]}:{x["rsid"]}') for x in leads["aad"] + leads[renal]]
        for i, (chrom, start, end, lead_ids) in enumerate(merge_windows(union), 1):
            window_id = f"{renal}_W{i:03d}"
            windows.append({"pair": f"aad__{renal}", "renal_trait": renal, "window_id": window_id, "chrom": chrom, "start_b37": start, "end_b37": end, "n_leads": len(lead_ids)})
            for lead in lead_ids:
                source_trait, rsid = lead.split(":", 1)
                memberships.append({"window_id": window_id, "source_trait": source_trait, "lead_rsid": rsid})

    for name, rows in [("coloc_windows.tsv", windows), ("window_lead_membership.tsv", memberships)]:
        with open(root / name, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    summary = {"lead_counts": {k: len(v) for k, v in leads.items()}, "window_counts": {renal: sum(w["renal_trait"] == renal for w in windows) for renal in RENAL}, "total_windows": len(windows)}
    with open(root / "locus_family_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
