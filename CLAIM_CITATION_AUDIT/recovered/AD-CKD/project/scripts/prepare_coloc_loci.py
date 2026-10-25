#!/usr/bin/env python3
"""Extract prespecified GWS rows and build merged lead-SNP windows."""

import csv
import gzip
import json
import math
from pathlib import Path


THRESHOLD = 5e-8


def normalize_rsid(value):
    value = (value or "").strip()
    return value if value.startswith("rs") and value[2:].isdigit() else None


def parse_p(value):
    try:
        p = float(value)
    except (TypeError, ValueError):
        return None
    return p if math.isfinite(p) and 0 <= p <= 1 else None


def merge_windows(rows, flank=500_000):
    intervals = sorted((str(chrom), max(1, int(pos) - flank), int(pos) + flank, lead) for chrom, pos, lead in rows)
    merged = []
    for chrom, start, end, lead in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            old = merged[-1]
            merged[-1] = (chrom, old[1], max(old[2], end), old[3] + [lead])
        else:
            merged.append((chrom, start, end, [lead]))
    return merged


SPECS = {
    "aad": dict(path="inputs/finngen_R13_I9_AORTDIS.gz", delimiter="\t", p="pval", rsid="rsids", chrom="#chrom", pos="pos", a1="alt", a2="ref", beta="beta", se="sebeta", freq="af_alt", n=None, build="GRCh38"),
    "egfrcrea": dict(path="inputs/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz", delimiter="\t", p="P.value.GC", rsid="RSID", chrom="chr", pos="pos", a1="Allele1", a2="Allele2", beta="Effect", se="StdErr.GC", freq="Freq1", n="n", build="GRCh37"),
    "egfrcys": dict(path="inputs/metal_eGFRcys_meta1.TBL.map.annot.gc.gz", delimiter="\t", p="P.value.GC", rsid="RSID", chrom="chr", pos="pos", a1="Allele1", a2="Allele2", beta="Effect", se="StdErr.GC", freq="Freq1", n="n", build="GRCh37"),
    "cystatin_c": dict(path="inputs/34226706-GCST90025945-EFO_0004617.h.tsv.gz", delimiter="\t", p="p_value", rsid="hm_rsid", chrom="hm_chrom", pos="hm_pos", a1="hm_effect_allele", a2="hm_other_allele", beta="hm_beta", se="standard_error", freq="hm_effect_allele_frequency", n=None, fixed_n="437846", build="GRCh37"),
    "bun": dict(path="inputs/BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt.gz", delimiter=None, p="P-value", rsid="RSID", chrom="Chr", pos="Pos_b37", a1="Allele1", a2="Allele2", beta="Effect", se="StdErr", freq="Freq1", n="n_total_sum", build="GRCh37"),
    "ckd": dict(path="inputs/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz", delimiter=None, p="P-value", rsid="RSID", chrom="Chr", pos="Pos_b37", a1="Allele1", a2="Allele2", beta="Effect", se="StdErr", freq="Freq1", n="n_total_sum", build="GRCh37"),
    "uacr": dict(path="inputs/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz", delimiter=None, p="P-value", rsid="RSID", chrom="Chr", pos="Pos_b37", a1="Allele1", a2="Allele2", beta="Effect", se="StdErr", freq="Freq1", n="n_total_sum", build="GRCh37"),
}


def split(line, delimiter):
    return line.rstrip("\n").split(delimiter) if delimiter else line.split()


def extract_trait(label, spec, output_dir):
    counts = {"rows": 0, "gws": 0, "written": 0, "invalid_rsid": 0}
    out = Path(output_dir) / f"{label}_gws.tsv"
    with gzip.open(spec["path"], "rt", errors="replace") as source, open(out, "w", newline="") as target:
        header = split(source.readline(), spec["delimiter"])
        idx = {name: header.index(name) for name in header}
        fields = ["trait", "rsid", "chrom", "pos", "build", "a1", "a2", "beta", "se", "p", "freq", "n"]
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            counts["rows"] += 1
            values = split(line, spec["delimiter"])
            p = parse_p(values[idx[spec["p"]]])
            if p is None or p >= THRESHOLD:
                continue
            counts["gws"] += 1
            rsid = normalize_rsid(values[idx[spec["rsid"]]])
            if rsid is None:
                counts["invalid_rsid"] += 1
                continue
            def get(key, default="NA"):
                col = spec.get(key)
                return values[idx[col]] if col else spec.get("fixed_" + key, default)
            writer.writerow(dict(trait=label, rsid=rsid, chrom=get("chrom"), pos=get("pos"), build=spec["build"], a1=get("a1").upper(), a2=get("a2").upper(), beta=get("beta"), se=get("se"), p=p, freq=get("freq"), n=get("n")))
            counts["written"] += 1
    with open(Path(output_dir) / f"{label}_clump_input.tsv", "w", newline="") as handle:
        handle.write("SNP\tP\n")
        with open(out) as source:
            next(source)
            for row in csv.DictReader(source, fieldnames=fields, delimiter="\t"):
                handle.write(f'{row["rsid"]}\t{row["p"]}\n')
    return counts


def main():
    output = Path("project/results/shared_genetics/coloc")
    output.mkdir(parents=True, exist_ok=True)
    summary = {}
    for label, spec in SPECS.items():
        summary[label] = extract_trait(label, spec, output)
        print(label, summary[label])
    with open(output / "gws_extraction_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
