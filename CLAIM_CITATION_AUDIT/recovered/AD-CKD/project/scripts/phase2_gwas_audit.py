#!/usr/bin/env python3
"""Audit Phase 2 GWAS inputs without loading whole summary-statistic files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


REQUIRED = {
    "GCST90103634_eGFR": [
        "metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",
        "GCST90103634.tsv.gz",
    ],
    "FinnGen_R13_I9_AORTDIS": ["finngen_R13_I9_AORTDIS.gz"],
}


def first_existing(directories: list[Path], names: list[str]) -> Path | None:
    for directory in directories:
        for name in names:
            candidate = directory / name
            if candidate.is_file():
                return candidate
    return None


def read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rt", errors="replace") as handle:
        return handle.readline().rstrip("\n\r").split("\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(
    path: Path | None,
    *,
    trait: str,
    build: str,
    ancestry: str,
    sample_size: str,
    required_fields: set[str],
    effect_allele: str,
    eligible: bool,
) -> dict:
    fields = read_header(path) if path else []
    missing_fields = sorted(required_fields - set(fields)) if path else sorted(required_fields)
    return {
        "path": str(path) if path else "MISSING",
        "present": bool(path),
        "trait": trait,
        "build": build,
        "ancestry": ancestry,
        "sample_size": sample_size,
        "fields": fields,
        "missing_fields": missing_fields,
        "field_status": "PASS" if path and not missing_fields else "FAIL",
        "effect_allele": effect_allele,
        "eligible": bool(path) and not missing_fields and eligible,
    }


def audit_inputs(gwas_dir: Path, downloads_dir: Path) -> dict:
    directories = [gwas_dir, downloads_dir]
    egfr = first_existing(directories, REQUIRED["GCST90103634_eGFR"])
    aad = first_existing(directories, REQUIRED["FinnGen_R13_I9_AORTDIS"])
    sbp = first_existing(directories, ["continuous-4080-both_sexes-irnt.tsv.bgz"])
    unknown_1 = first_existing(directories, ["GCST90469223.tsv.gz"])
    unknown_2 = first_existing(directories, ["GCST90478526.tsv.gz"])
    ld_prefix = gwas_dir / "ld" / "EUR"
    ld_present = all(Path(f"{ld_prefix}.{ext}").is_file() for ext in ("bed", "bim", "fam"))

    inputs = {
        "GCST90103634_eGFR": record(
            egfr,
            trait="eGFRcrea (Stanzick 2021)",
            build="GRCh37",
            ancestry="European",
            sample_size="1004040",
            required_fields={"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "RSID"},
            effect_allele="Allele1/coded allele; confirm against file README",
            eligible=True,
        ),
        "FinnGen_R13_I9_AORTDIS": record(
            aad,
            trait="Aortic dissection (I9_AORTDIS)",
            build="GRCh38",
            ancestry="Finnish/European",
            sample_size="capture from R13 endpoint metadata",
            required_fields={"#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta"},
            effect_allele="alt",
            eligible=True,
        ),
        "UKB_4080_SBP": record(
            sbp,
            trait="Systolic blood pressure, automated reading, IRNT",
            build="GRCh37 (Pan-UKB flat files; recheck before harmonisation)",
            ancestry="multi-ancestry file; use EUR-specific columns only",
            sample_size="416955 full cohort; 396663 EUR",
            required_fields={"chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR"},
            effect_allele="alt",
            eligible=True,
        ),
        "GCST90469223": record(
            unknown_1,
            trait="UNKNOWN",
            build="UNKNOWN",
            ancestry="UNKNOWN",
            sample_size="47745 observed in file; provenance unverified",
            required_fields=set(),
            effect_allele="effect_allele according to header only",
            eligible=False,
        ),
        "GCST90478526": record(
            unknown_2,
            trait="UNKNOWN",
            build="UNKNOWN",
            ancestry="UNKNOWN",
            sample_size="3469 cases + 446868 controls observed; provenance unverified",
            required_fields=set(),
            effect_allele="effect_allele according to header only",
            eligible=False,
        ),
        "1000G_v3_EUR": {
            "path": str(ld_prefix),
            "present": ld_present,
            "trait": "LD reference panel",
            "build": "GRCh37",
            "ancestry": "European",
            "sample_size": "503 individuals (verify from FAM)",
            "fields": ["bed", "bim", "fam"],
            "missing_fields": [] if ld_present else ["bed/bim/fam"],
            "field_status": "PASS" if ld_present else "FAIL",
            "effect_allele": "PLINK BIM A1/A2; not an MR effect estimate",
            "eligible": ld_present,
        },
    }
    missing = [name for name in REQUIRED if not inputs[name]["present"]]
    return {
        "schema": "ad-ckd.phase2.gwas-audit.v1",
        "mr_ready": not missing and all(inputs[name]["eligible"] for name in REQUIRED),
        "missing_required": missing,
        "inputs": inputs,
        "interpretation_guardrail": "MVMR estimates conditional direct effects; it does not by itself prove mediation.",
    }


def write_outputs(result: dict, out_dir: Path, checksum: bool) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    if checksum:
        for item in result["inputs"].values():
            path = Path(item["path"])
            if item["present"] and path.is_file():
                item["sha256"] = sha256(path)

    (out_dir / "phase2_gwas_audit.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    columns = [
        "dataset", "present", "eligible", "trait", "build", "ancestry",
        "sample_size", "field_status", "effect_allele", "path", "sha256",
    ]
    with (out_dir / "mr_input_registry.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for name, item in result["inputs"].items():
            writer.writerow({"dataset": name, **{key: item.get(key, "") for key in columns[1:]}})

    state = "READY" if result["mr_ready"] else "BLOCKED_MISSING_PRIMARY_GWAS"
    missing = ", ".join(result["missing_required"]) or "none"
    lines = [
        "# Phase 2 MR 输入审计", "", f"状态：**{state}**", "",
        f"缺失的必需主文件：{missing}", "",
        "当前冻结规则：只有 trait、build、ancestry、样本量、必要字段和效应等位基因方向均可验证的文件才可进入 MR。",
        "未知 trait 的 GCST 文件保持隔离，不得因字段相似而替代 eGFR 或主动脉夹层数据。",
        "SBP 仅使用 EUR 专属的 `beta_EUR`、`se_EUR`、`af_EUR`，效应等位基因为 `alt`。",
        "MVMR 只估计条件直接效应，不能单独证明中介。", "",
        "机器可读结果：`phase2_gwas_audit.json`；冻结表：`mr_input_registry.tsv`。", "",
    ]
    (out_dir / "PHASE2_INPUT_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-dir", type=Path, required=True)
    parser.add_argument("--downloads-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--checksum", action="store_true")
    args = parser.parse_args()
    result = audit_inputs(args.gwas_dir, args.downloads_dir)
    write_outputs(result, args.out_dir, args.checksum)
    print(json.dumps({"mr_ready": result["mr_ready"], "missing_required": result["missing_required"]}))
    return 0 if result["mr_ready"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
