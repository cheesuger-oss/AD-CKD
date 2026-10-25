#!/usr/bin/env python3
"""Extract and harmonize all variants in the frozen colocalization windows."""

import bisect
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

from prepare_coloc_loci import SPECS, normalize_rsid, parse_p, split


def parse_values(line, delimiter, expected_columns):
    values = split(line, delimiter)
    return values if len(values) == expected_columns else None


class WindowIndex:
    def __init__(self, windows):
        self.by_chrom = defaultdict(list)
        for chrom, start, end, window_id in windows:
            self.by_chrom[str(chrom)].append((int(start), int(end), window_id))
        self.starts = {}
        for chrom in self.by_chrom:
            self.by_chrom[chrom].sort()
            self.starts[chrom] = [x[0] for x in self.by_chrom[chrom]]

    def find(self, chrom, pos):
        chrom = str(chrom); pos = int(pos)
        if chrom not in self.by_chrom:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos) - 1
        if i >= 0:
            start, end, window_id = self.by_chrom[chrom][i]
            if start <= pos <= end:
                return window_id
        return None


def finite(value, positive=False):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x) or (positive and x <= 0):
        return None
    return x


def harmonize(renal, aad):
    if {renal["a1"], renal["a2"]} != {aad["a1"], aad["a2"]}:
        return None
    rb, ab = finite(renal["beta"]), finite(aad["beta"])
    rse, ase = finite(renal["se"], True), finite(aad["se"], True)
    rp, ap = parse_p(renal["p"]), parse_p(aad["p"])
    rf, af = finite(renal["freq"]), finite(aad["freq"])
    rn = finite(renal["n"], True)
    if None in (rb, ab, rse, ase, rp, ap, rf, af, rn) or not (0 < rf < 1 and 0 < af < 1):
        return None
    if aad["a1"] != renal["a1"]:
        ab, af = -ab, 1 - af
    return {"rsid": renal["rsid"], "a1": renal["a1"], "a2": renal["a2"], "renal_beta": rb,
            "renal_varbeta": rse * rse, "renal_p": rp, "renal_freq": rf, "renal_n": rn,
            "aad_beta": ab, "aad_varbeta": ase * ase, "aad_p": ap, "aad_freq": af,
            "aad_n": 519870.0, "aad_case_fraction": 1262 / 519870}


def read_windows(root):
    grouped = defaultdict(list)
    with open(root / "coloc_windows.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["renal_trait"]].append((row["chrom"], int(row["start_b37"]), int(row["end_b37"]), row["window_id"]))
    return grouped


def standard_row(values, idx, spec):
    def get(key, default="NA"):
        col = spec.get(key)
        return values[idx[col]] if col else spec.get("fixed_" + key, default)
    return {"rsid": normalize_rsid(get("rsid")), "chrom": get("chrom"), "pos": get("pos"),
            "a1": get("a1").upper(), "a2": get("a2").upper(), "beta": get("beta"), "se": get("se"),
            "p": get("p"), "freq": get("freq"), "n": get("n")}


def extract_renal(trait, spec, index, root, wanted):
    counts = {"rows": 0, "malformed_rows": 0, "in_windows": 0, "with_rsid": 0}
    out = root / f"{trait}_regional.tsv.gz"
    fields = ["window_id", "rsid", "chrom", "pos", "a1", "a2", "beta", "se", "p", "freq", "n"]
    with gzip.open(spec["path"], "rt", errors="replace") as source, gzip.open(out, "wt", newline="") as target:
        header = split(source.readline(), spec["delimiter"]); idx = {x: header.index(x) for x in header}
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for line in source:
            counts["rows"] += 1
            values = parse_values(line, spec["delimiter"], len(header))
            if values is None:
                counts["malformed_rows"] += 1
                continue
            row = standard_row(values, idx, spec)
            window_id = index.find(row["chrom"], row["pos"])
            if not window_id: continue
            counts["in_windows"] += 1
            if not row["rsid"]: continue
            counts["with_rsid"] += 1; row["window_id"] = window_id; writer.writerow(row)
            wanted[row["rsid"]].append((trait, window_id))
    return counts


def extract_aad(spec, root, wanted):
    fields = ["trait", "window_id", "rsid", "a1", "a2", "beta", "se", "p", "freq"]
    out = root / "aad_regional_by_pair.tsv.gz"; counts = {"rows": 0, "malformed_rows": 0, "matched_rows": 0, "written_memberships": 0}
    with gzip.open(spec["path"], "rt", errors="replace") as source, gzip.open(out, "wt", newline="") as target:
        header = split(source.readline(), spec["delimiter"]); idx = {x: header.index(x) for x in header}
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for line in source:
            counts["rows"] += 1; values = parse_values(line, spec["delimiter"], len(header))
            if values is None:
                counts["malformed_rows"] += 1
                continue
            row = standard_row(values, idx, spec)
            memberships = wanted.get(row["rsid"])
            if not memberships: continue
            counts["matched_rows"] += 1
            for trait, window_id in memberships:
                writer.writerow({"trait": trait, "window_id": window_id, "rsid": row["rsid"], "a1": row["a1"], "a2": row["a2"], "beta": row["beta"], "se": row["se"], "p": row["p"], "freq": row["freq"]})
                counts["written_memberships"] += 1
    return counts


def join_all(root, traits):
    aad = defaultdict(dict)
    with gzip.open(root / "aad_regional_by_pair.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            aad[row["trait"]][(row["window_id"], row["rsid"])] = row
    summary = {}
    fields = ["window_id", "rsid", "a1", "a2", "renal_beta", "renal_varbeta", "renal_p", "renal_freq", "renal_n", "aad_beta", "aad_varbeta", "aad_p", "aad_freq", "aad_n", "aad_case_fraction"]
    for trait in traits:
        stats = defaultdict(lambda: {"renal_rows": 0, "aad_overlap": 0, "harmonized": 0})
        with gzip.open(root / f"{trait}_regional.tsv.gz", "rt") as source, gzip.open(root / f"aad__{trait}_harmonized.tsv.gz", "wt", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
            for renal in csv.DictReader(source, delimiter="\t"):
                st = stats[renal["window_id"]]; st["renal_rows"] += 1
                arow = aad[trait].get((renal["window_id"], renal["rsid"]))
                if not arow: continue
                st["aad_overlap"] += 1; joined = harmonize(renal, arow)
                if joined is None: continue
                st["harmonized"] += 1; joined["window_id"] = renal["window_id"]; writer.writerow(joined)
        summary[trait] = dict(stats)
    with open(root / "regional_identifiability.json", "w") as handle: json.dump(summary, handle, indent=2)
    return summary


def main():
    root = Path("project/results/shared_genetics/coloc"); windows = read_windows(root); wanted = defaultdict(list); summary = {"renal": {}}
    for trait, regions in windows.items():
        summary["renal"][trait] = extract_renal(trait, SPECS[trait], WindowIndex(regions), root, wanted)
        print(trait, summary["renal"][trait], flush=True)
    summary["aad"] = extract_aad(SPECS["aad"], root, wanted); print("aad", summary["aad"], flush=True)
    ident = join_all(root, list(windows)); summary["identifiable_windows_ge100"] = {t: sum(v["harmonized"] >= 100 for v in ident[t].values()) for t in ident}
    with open(root / "regional_extraction_summary.json", "w") as handle: json.dump(summary, handle, indent=2)
    print(summary["identifiable_windows_ge100"])


if __name__ == "__main__":
    main()
