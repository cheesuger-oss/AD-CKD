#!/usr/bin/env python3
"""Prepare allele-checked Pan-UKB SBP candidates for joint MVMR clumping."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _alleles_match(ref: str, alt: str, a1: str, a2: str) -> bool:
    observed = {ref.upper(), alt.upper()}
    direct = {a1.upper(), a2.upper()}
    complement = {a1.upper().translate(COMPLEMENT), a2.upper().translate(COMPLEMENT)}
    return observed == direct or observed == complement


def map_sbp_to_reference(sbp_file: Path, bim_file: Path, output: Path) -> dict:
    with sbp_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("SBP candidate file has no header")
        fieldnames = reader.fieldnames
        sbp_rows = list(reader)
    wanted_positions = {(row["chr"], row["pos"]) for row in sbp_rows}
    reference = defaultdict(list)
    with bim_file.open() as handle:
        for line in handle:
            chrom, rsid, _, pos, a1, a2 = line.rstrip("\n").split()[:6]
            if (chrom, pos) in wanted_positions:
                reference[(chrom, pos)].append((rsid, a1, a2))

    input_rows = mapped_rows = ambiguous_rows = 0
    with output.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=fieldnames + ["RSID"], delimiter="\t")
        writer.writeheader()
        for row in sbp_rows:
            input_rows += 1
            matches = {
                rsid
                for rsid, a1, a2 in reference.get((row["chr"], row["pos"]), [])
                if _alleles_match(row["ref"], row["alt"], a1, a2)
            }
            if len(matches) == 1:
                row["RSID"] = next(iter(matches))
                writer.writerow(row)
                mapped_rows += 1
            elif len(matches) > 1:
                ambiguous_rows += 1
    return {
        "input_rows": input_rows,
        "mapped_rows": mapped_rows,
        "ambiguous_rows": ambiguous_rows,
    }


def build_joint_clump_input(egfr_file: Path, sbp_file: Path, output: Path) -> dict:
    variants = {}
    egfr_candidates = sbp_candidates = 0
    with egfr_file.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            egfr_candidates += 1
            p_value = float(row["P.value.GC"])
            current = variants.get(row["RSID"])
            if current is None or p_value < current[0]:
                variants[row["RSID"]] = (p_value, row["chr"], row["pos"])
    with sbp_file.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sbp_candidates += 1
            p_value = 10 ** (-float(row["neglog10_pval_EUR"]))
            if not math.isfinite(p_value):
                continue
            current = variants.get(row["RSID"])
            if current is None or p_value < current[0]:
                variants[row["RSID"]] = (p_value, row["chr"], row["pos"])
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["SNP", "P", "CHR", "BP"])
        for rsid, (p_value, chrom, pos) in sorted(variants.items(), key=lambda item: item[1][0]):
            writer.writerow([rsid, f"{p_value:.12g}", chrom, pos])
    return {
        "egfr_candidates": egfr_candidates,
        "sbp_candidates": sbp_candidates,
        "union_variants": len(variants),
    }


def export_clumped_targets(
    clumped_file: Path, joint_file: Path, ids_output: Path, positions_output: Path
) -> int:
    selected = []
    with clumped_file.open() as handle:
        header = handle.readline().split()
        snp_index = header.index("SNP")
        for line in handle:
            fields = line.split()
            if fields:
                selected.append(fields[snp_index])
    coordinates = {}
    with joint_file.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            coordinates[row["SNP"]] = (row["CHR"], row["BP"])
    missing = [rsid for rsid in selected if rsid not in coordinates]
    if missing:
        raise ValueError(f"Clumped variants absent from joint input: {missing[:3]}")
    ids_output.write_text("".join(f"{rsid}\n" for rsid in selected))
    with positions_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chr", "pos"])
        writer.writerows(coordinates[rsid] for rsid in selected)
    return len(selected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint", action="store_true")
    parser.add_argument("--export-clumped", action="store_true")
    parser.add_argument("sbp_file", type=Path, nargs="?")
    parser.add_argument("bim_file", type=Path, nargs="?")
    parser.add_argument("output", type=Path, nargs="?")
    parser.add_argument("--egfr-file", type=Path)
    parser.add_argument("--clumped-file", type=Path)
    parser.add_argument("--ids-output", type=Path)
    parser.add_argument("--positions-output", type=Path)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    if args.export_clumped:
        required = [args.clumped_file, args.sbp_file, args.ids_output, args.positions_output]
        if any(value is None for value in required):
            parser.error("--export-clumped requires joint input positional plus --clumped-file, --ids-output and --positions-output")
        summary = {"clumped_instruments": export_clumped_targets(
            args.clumped_file, args.sbp_file, args.ids_output, args.positions_output
        )}
    elif args.joint:
        if args.egfr_file is None:
            parser.error("--joint requires --egfr-file")
        summary = build_joint_clump_input(args.egfr_file, args.sbp_file, args.output)
    else:
        if any(value is None for value in [args.sbp_file, args.bim_file, args.output]):
            parser.error("mapping requires sbp_file, bim_file and output")
        summary = map_sbp_to_reference(args.sbp_file, args.bim_file, args.output)
    if args.summary:
        args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
