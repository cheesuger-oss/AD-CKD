#!/usr/bin/env python3

import csv
import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "project/scripts/phase3_freeze_union_genes.py"
OUT = ROOT / "project/results/network"


subprocess.run([sys.executable, str(SCRIPT), "--root", str(ROOT)], check=True)
union_path = OUT / "union_genes.txt"
genes = union_path.read_text(encoding="utf-8").splitlines()
assert genes == sorted(set(genes))
assert len(genes) == 579

with (OUT / "gene_source_audit.tsv").open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(rows) == len(genes)
assert list(rows[0]) == ["gene", "phase1_intersection", "omim_protocol_seed", "uremic_protocol_seed"]
assert sum(int(row["phase1_intersection"]) for row in rows) == 549
assert sum(int(row["omim_protocol_seed"]) for row in rows) == 20
assert sum(int(row["uremic_protocol_seed"]) for row in rows) == 26

expected = hashlib.sha256(union_path.read_bytes()).hexdigest()
recorded = (OUT / "union_genes.sha256").read_text(encoding="utf-8").split()[0]
assert recorded == expected
print("phase3 union freeze tests: PASS")
