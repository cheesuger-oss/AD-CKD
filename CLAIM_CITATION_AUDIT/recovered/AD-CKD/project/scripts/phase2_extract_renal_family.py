#!/usr/bin/env python3
"""Stream registered renal GWAS files into standardized genome-wide-significant candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


OUT_FIELDS = [
    "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval",
    "n", "chr", "pos", "exposure_id", "trait", "scale",
]


def split_line(line: str, delimiter: str) -> list[str]:
    return line.rstrip("\r\n").split("\t") if delimiter == "tab" else line.split()


def extract(row: dict[str, str], output: Path, p_threshold: float = 5e-8) -> dict[str, object]:
    source = Path(row["source_file"])
    if row["status"] != "VERIFIED":
        return {"exposure_id": row["exposure_id"], "status": row["status"], "source": str(source)}
    best: dict[str, dict[str, object]] = {}
    scanned = valid = 0
    with gzip.open(source, "rt", encoding="utf-8", errors="replace") as handle:
        columns = split_line(handle.readline(), row["delimiter"])
        field_keys = ("chr_col", "pos_col", "snp_col", "ea_col", "oa_col", "eaf_col",
                      "effect_col", "se_col", "p_col")
        required = [row[name] for name in field_keys]
        fixed_n = None
        try:
            fixed_n = float(row["n_col"])
        except ValueError:
            required.append(row["n_col"])
        missing = [field for field in required if field not in columns]
        if missing:
            raise ValueError(f"{row['exposure_id']}: missing fields {missing}; columns={columns}")
        idx = {field: columns.index(field) for field in required}
        for line in handle:
            scanned += 1
            values = split_line(line, row["delimiter"])
            if len(values) < len(columns):
                continue
            try:
                snp = values[idx[row["snp_col"]]]
                pval = float(values[idx[row["p_col"]]])
                raw_effect = float(values[idx[row["effect_col"]]])
                se = float(values[idx[row["se_col"]]])
                eaf = float(values[idx[row["eaf_col"]]])
                n = fixed_n if fixed_n is not None else float(values[idx[row["n_col"]]])
                chrom = values[idx[row["chr_col"]]]
                pos = values[idx[row["pos_col"]]]
                ea = values[idx[row["ea_col"]]].upper()
                oa = values[idx[row["oa_col"]]].upper()
            except (ValueError, IndexError):
                continue
            if not snp.startswith("rs") or not (0 < pval < p_threshold) or se <= 0 or not (0 <= eaf <= 1):
                continue
            beta = math.log(raw_effect) if row["scale"] == "odds_ratio" else raw_effect
            if not all(math.isfinite(x) for x in (beta, se, pval, eaf, n)):
                continue
            record = {
                "SNP": snp, "beta": beta, "se": se, "effect_allele": ea, "other_allele": oa,
                "eaf": eaf, "pval": pval, "n": n, "chr": chrom, "pos": pos,
                "exposure_id": row["exposure_id"], "trait": row["trait"], "scale": row["scale"],
            }
            valid += 1
            if snp not in best or pval < float(best[snp]["pval"]):
                best[snp] = record
    records = sorted(best.values(), key=lambda item: (float(item["pval"]), str(item["SNP"])))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    return {
        "exposure_id": row["exposure_id"], "status": "COMPLETE", "source": str(source),
        "rows_scanned": scanned, "significant_rows": valid, "unique_significant_rsids": len(records),
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--exposure")
    args = parser.parse_args()
    root = args.root.resolve()
    registry = root / "project/results/mr/renal_exposure_registry.tsv"
    outdir = root / "project/data/gwas/processed/renal_family"
    summaries = []
    with registry.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if args.exposure and row["exposure_id"] != args.exposure:
                continue
            print(f"scanning {row['exposure_id']}", flush=True)
            summary = extract(row, outdir / f"{row['exposure_id']}_gws_candidates.tsv")
            summaries.append(summary)
            print(json.dumps(summary, ensure_ascii=False), flush=True)
    (outdir / "extraction_summary.json").write_text(
        json.dumps(summaries, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
