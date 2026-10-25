#!/usr/bin/env python3
"""Repair serum-cystatin-C regional membership using reference rsID positions."""

import csv
import gzip
from collections import defaultdict
from pathlib import Path

from extract_coloc_regions import WindowIndex, standard_row, harmonize, parse_values
from prepare_coloc_loci import SPECS, split


ROOT = Path("project/results/shared_genetics/coloc")


def load_cystatin_reference_windows():
    windows=[]
    with open(ROOT/"coloc_windows.tsv") as h:
        for r in csv.DictReader(h,delimiter="\t"):
            if r["renal_trait"]=="cystatin_c": windows.append((r["chrom"],r["start_b37"],r["end_b37"],r["window_id"]))
    index=WindowIndex(windows); mapping={}
    with open("project/data/gwas/ld/EUR.bim") as h:
        for line in h:
            chrom,rsid,_,pos,*_=line.split()
            wid=index.find(chrom,pos)
            if wid: mapping[rsid]=wid
    return mapping


def main():
    mapping=load_cystatin_reference_windows(); spec=SPECS["cystatin_c"]
    renal_out=ROOT/"cystatin_c_regional_b37rsid.tsv.gz"; wanted=defaultdict(list)
    fields=["window_id","rsid","chrom","pos","a1","a2","beta","se","p","freq","n"]
    counts={"reference_rsids_in_windows":len(mapping),"source_rows":0,"mapped_rows":0}
    with gzip.open(spec["path"],"rt",errors="replace") as source,gzip.open(renal_out,"wt",newline="") as target:
        header=split(source.readline(),spec["delimiter"]); idx={x:header.index(x) for x in header}
        writer=csv.DictWriter(target,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader()
        for line in source:
            counts["source_rows"]+=1; values=parse_values(line,spec["delimiter"],len(header))
            if values is None: continue
            row=standard_row(values,idx,spec); wid=mapping.get(row["rsid"])
            if not wid: continue
            row["window_id"]=wid; row["chrom"]="NA"; row["pos"]="NA"; writer.writerow(row)
            wanted[row["rsid"]].append(wid); counts["mapped_rows"]+=1

    aad_spec=SPECS["aad"]; aad=defaultdict(dict); counts.update(aad_rows=0,aad_matched_rows=0)
    with gzip.open(aad_spec["path"],"rt",errors="replace") as source:
        header=split(source.readline(),aad_spec["delimiter"]); idx={x:header.index(x) for x in header}
        for line in source:
            counts["aad_rows"]+=1; values=parse_values(line,aad_spec["delimiter"],len(header))
            if values is None: continue
            row=standard_row(values,idx,aad_spec); wids=wanted.get(row["rsid"])
            if not wids: continue
            counts["aad_matched_rows"]+=1
            for wid in wids: aad[(wid,row["rsid"])]=row

    out=ROOT/"aad__cystatin_c_harmonized_b37rsid.tsv.gz"
    outfields=["window_id","rsid","a1","a2","renal_beta","renal_varbeta","renal_p","renal_freq","renal_n","aad_beta","aad_varbeta","aad_p","aad_freq","aad_n","aad_case_fraction"]
    stats=defaultdict(lambda:{"renal_rows":0,"aad_overlap":0,"harmonized":0})
    with gzip.open(renal_out,"rt") as source,gzip.open(out,"wt",newline="") as target:
        writer=csv.DictWriter(target,fieldnames=outfields,delimiter="\t",lineterminator="\n");writer.writeheader()
        for renal in csv.DictReader(source,delimiter="\t"):
            s=stats[renal["window_id"]];s["renal_rows"]+=1; ar=aad.get((renal["window_id"],renal["rsid"]))
            if not ar: continue
            s["aad_overlap"]+=1; joined=harmonize(renal,ar)
            if joined is None: continue
            s["harmonized"]+=1;joined["window_id"]=renal["window_id"];writer.writerow(joined)
    counts["windows_ge100_harmonized"]=sum(s["harmonized"]>=100 for s in stats.values())
    counts["windows_with_rows"]=len(stats)
    print(counts)


if __name__=="__main__": main()
