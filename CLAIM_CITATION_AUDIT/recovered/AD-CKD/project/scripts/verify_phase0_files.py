#!/usr/bin/env python3
import gzip
import hashlib
import json
from pathlib import Path
import sys
import zipfile

root = Path(__file__).resolve().parents[2]
failures = []
checks = 0


def check(condition, label):
    global checks
    checks += 1
    print(("PASS" if condition else "FAIL") + "\t" + label)
    if not condition:
        failures.append(label)


for relative in (
    "project/data/geo",
    "project/data/gwas",
    "project/data/string",
    "project/scripts",
    "project/results/geo",
    "project/results/mr",
    "project/results/network",
    "project/results/virtual_perturbation",
    "project/results/cohort",
    "project/manuscript",
):
    check((root / relative).is_dir(), f"directory {relative}")

geo_files = {
    "project/data/geo/GSE146638_Effect_of_CKD_on_calcified_aorta_by_RNAseq_Processed_data.xlsx":
        "eab101f6c0d0e460d4ff198fd3055f29f15bdfb6e5004d7afc12501ba93cb108",
    "project/data/geo/GSE325608_raw_counts.txt.gz":
        "a0016ec31daff4865d33b206ecedd513f871654f5699310b66e47b2cc228ec1a",
}
for relative, expected in geo_files.items():
    path = root / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else ""
    check(actual == expected, f"checksum {relative}")

xlsx = root / next(iter(geo_files))
try:
    with zipfile.ZipFile(xlsx) as archive:
        check(archive.testzip() is None, "GSE146638 XLSX archive integrity")
except Exception:
    check(False, "GSE146638 XLSX archive integrity")

try:
    with gzip.open(root / "project/data/geo/GSE325608_raw_counts.txt.gz", "rb") as handle:
        while handle.read(1024 * 1024):
            pass
    check(True, "GSE325608 gzip integrity")
except Exception:
    check(False, "GSE325608 gzip integrity")

sources = (root / "project/data/data_sources.tsv").read_text(encoding="utf-8")
for token in ("GSE146638", "GSE325608", "GCST90103634", "FinnGen_R13_I9_AORTDIS", "STRING_v12", "UKB_4080_IRNT"):
    check(token in sources, f"registered source {token}")

lock_path = root / "renv.lock"
check(lock_path.is_file(), "renv.lock exists")
lock = json.loads(lock_path.read_text(encoding="utf-8")) if lock_path.is_file() else {"Packages": {}}
required = (
    "renv", "GEOquery", "DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db",
    "biomaRt", "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "ggplot2",
    "pheatmap", "stringr", "dplyr", "data.table", "igraph",
)
library = root / "renv/library/macos/R-4.5/aarch64-apple-darwin20"
bootstrap = root / "project/Rlib-4.5"
for package in required:
    description = (library / package / "DESCRIPTION")
    if not description.is_file():
        description = bootstrap / package / "DESCRIPTION"
    version = ""
    if description.is_file():
        for line in description.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("Version: "):
                version = line.split(": ", 1)[1]
                break
    locked = lock.get("Packages", {}).get(package, {}).get("Version", "")
    check(bool(version) and version == locked, f"lock matches {package} {version or 'missing'}")

check((root / "project/venv/bin/python").is_file(), "Python virtual environment")
check((root / ".light/project_card.md").is_file(), "project memory card")
check((root / ".light/decision_log.md").is_file(), "decision log")

print(f"SUMMARY\tpasses_or_checks={checks - len(failures)}/{checks}\tfailures={len(failures)}")
sys.exit(1 if failures else 0)
