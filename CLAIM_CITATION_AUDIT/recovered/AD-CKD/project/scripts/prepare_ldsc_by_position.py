#!/usr/bin/env python3
"""Map a position-only GWAS file to HapMap3 rsIDs for LDSC munging."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def choose_match(candidates, effect, other):
    matches = []
    for rsid, a1, a2 in candidates:
        if {effect, other} == {a1, a2}:
            matches.append((rsid, effect, other))
    return matches[0] if len(matches) == 1 else None


def load_reference(ld_dir):
    alleles = {}
    with open(Path(ld_dir) / "w_hm3.snplist", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            alleles[row["SNP"]] = (row["A1"], row["A2"])
    by_pos = defaultdict(list)
    for chrom in range(1, 23):
        with gzip.open(Path(ld_dir) / f"{chrom}.l2.ldscore.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rsid = row["SNP"]
                if rsid in alleles:
                    by_pos[(row["CHR"], row["BP"])].append((rsid, *alleles[rsid]))
    return by_pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--ld-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-size", required=True, type=float)
    args = parser.parse_args()

    reference = load_reference(args.ld_dir)
    counts = {"read": 0, "position_match": 0, "written": 0, "ambiguous_or_allele_mismatch": 0}
    opener = gzip.open if args.input.endswith(".gz") else open
    with opener(args.input, "rt", newline="") as source, gzip.open(args.output, "wt", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        fields = ["SNP", "A1", "A2", "P", "BETA", "N"]
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            counts["read"] += 1
            candidates = reference.get((row["chromosome"], row["base_pair_location"]))
            if not candidates:
                continue
            counts["position_match"] += 1
            match = choose_match(candidates, row["effect_allele"], row["other_allele"])
            if match is None:
                counts["ambiguous_or_allele_mismatch"] += 1
                continue
            rsid, a1, a2 = match
            writer.writerow({"SNP": rsid, "A1": a1, "A2": a2, "P": row["p_value"], "BETA": row["beta"], "N": args.sample_size})
            counts["written"] += 1
    print(counts)


if __name__ == "__main__":
    main()
