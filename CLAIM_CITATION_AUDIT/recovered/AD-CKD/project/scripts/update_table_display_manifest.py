#!/usr/bin/env python3
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "project/manuscript/DISPLAY_ITEM_MANIFEST.tsv"

def sha(rel):
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()

def joined_hashes(paths):
    return ";".join(f"{p}={sha(p)}" for p in paths)

with MANIFEST.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)

main_book = "project/manuscript/outputs/ad_ckd_submission/Main_Tables_1-2.xlsx"
main_gen = "project/scripts/build_main_tables.mjs"
supp_book = "project/manuscript/outputs/ad_ckd_submission/Supplementary_Tables_S1-S4.xlsx"
supp_gen = "project/scripts/build_supplementary_tables.mjs"

table1_sources = [
    "project/results/geo/GSE146638_sample_metadata.tsv",
    "project/results/geo/GSE325608_sample_metadata.tsv",
    "project/results/geo/GSE146638_DEG.csv",
    "project/results/geo/GSE325608_DEG.csv",
    "project/results/geo/GSE146638_human_1to1_DEG.csv",
    "project/results/geo/GSE325608_human_1to1_DEG.csv",
    "project/manuscript/METHODS_MAIN_APPROVED_SECTIONS.md",
]

for row in rows:
    item = row["item_id"]
    if item in {"Table 1", "Table 2"}:
        row["status"] = "FROZEN_PASS_WORKBOOK"
        row["visual_files"] = main_book
        row["visual_sha256"] = f"{main_book}={sha(main_book)}"
        row["generator_files"] = main_gen
        row["generator_sha256"] = f"{main_gen}={sha(main_gen)}"
        if item == "Table 1":
            row["source_data_files"] = ";".join(table1_sources)
            row["source_data_sha256"] = joined_hashes(table1_sources)
    elif item in {f"Supplementary Table S{i}" for i in range(1, 5)}:
        row["status"] = "FROZEN_PASS_WORKBOOK"
        row["visual_files"] = supp_book
        row["visual_sha256"] = f"{supp_book}={sha(supp_book)}"
        row["generator_files"] = supp_gen
        row["generator_sha256"] = f"{supp_gen}={sha(supp_gen)}"

with MANIFEST.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_ALL, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
