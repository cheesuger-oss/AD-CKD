#!/usr/bin/env python3
"""Summarize the frozen coloc family and render dependency-free SVG figures."""

import csv
import html
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path("project/results/shared_genetics/coloc")


def classify(row):
    h3, h4, ratio = (float(row[x]) for x in ("PP.H3", "PP.H4", "H4_over_H3H4"))
    if h4 >= 0.8 and ratio >= 0.8:
        return "PRIMARY_THRESHOLD_PRIOR_SENSITIVE"
    if h3 >= 0.8:
        return "DISTINCT_VARIANTS"
    return "NO_PRIMARY_SUPPORT"


def render_overview(rows, path):
    traits = ["egfrcrea", "egfrcys", "cystatin_c", "bun", "ckd", "uacr"]
    width, height, left, right, top = 1050, 520, 105, 35, 45
    panel_w = (width-left-right)/len(traits); plot_h = 360
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><rect width="100%" height="100%" fill="white"/>',
           '<style>text{font-family:Arial,sans-serif;fill:#222}.lab{font-size:15px}.small{font-size:12px}</style>']
    for j,t in enumerate(traits):
        xrows=[r for r in rows if r["trait"]==t]; x0=left+j*panel_w
        for i,r in enumerate(xrows):
            x=x0+(i+0.5)/len(xrows)*panel_w; y=top+(1-float(r["PP.H4"]))*plot_h
            status=classify(r); color="#D55E00" if status.startswith("PRIMARY") else ("#0072B2" if status=="DISTINCT_VARIANTS" else "#999999")
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="{color}" fill-opacity="0.75"/>')
        parts.append(f'<text class="small" x="{x0+panel_w/2:.1f}" y="{top+plot_h+24}" text-anchor="middle">{html.escape(t)}</text>')
        if j: parts.append(f'<line x1="{x0:.1f}" y1="{top}" x2="{x0:.1f}" y2="{top+plot_h}" stroke="#eee"/>')
    y80=top+0.2*plot_h; parts.append(f'<line x1="{left}" y1="{y80}" x2="{width-right}" y2="{y80}" stroke="#D55E00" stroke-dasharray="5,4"/>')
    for v in [0,0.2,0.4,0.6,0.8,1.0]:
        y=top+(1-v)*plot_h; parts.append(f'<text class="small" x="{left-12}" y="{y+4:.1f}" text-anchor="end">{v:.1f}</text>')
    parts += [f'<text class="lab" x="25" y="{top+plot_h/2}" transform="rotate(-90 25 {top+plot_h/2})" text-anchor="middle">PP.H4 (primary prior)</text>',
              f'<text class="lab" x="{left+(width-left-right)/2}" y="{height-28}" text-anchor="middle">Prespecified windows within each renal trait</text>', '</svg>']
    path.write_text("\n".join(parts)+"\n")


def render_sensitivity(all_rows, candidates, path):
    width,height,left,top=760,430,90,45; plot_w,plot_h=610,280
    colors={c:"#D55E00" if "cystatin" in c else "#0072B2" for c in candidates}
    pvals=[1e-6,5e-6,1e-5,5e-5]; xpos={p:left+i*plot_w/(len(pvals)-1) for i,p in enumerate(pvals)}
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><rect width="100%" height="100%" fill="white"/>','<style>text{font-family:Arial,sans-serif;fill:#222}.lab{font-size:15px}.small{font-size:12px}</style>']
    for cand in candidates:
        vals=sorted([r for r in all_rows if f'{r["trait"]}:{r["window_id"]}'==cand],key=lambda r:float(r['p12']))
        pts=[(xpos[float(r['p12'])],top+(1-float(r['PP.H4']))*plot_h) for r in vals]
        parts.append('<polyline points="'+' '.join(f'{x:.1f},{y:.1f}' for x,y in pts)+f'" fill="none" stroke="{colors[cand]}" stroke-width="3"/>')
        for x,y in pts: parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{colors[cand]}"/>')
    y80=top+.2*plot_h; parts.append(f'<line x1="{left}" y1="{y80}" x2="{left+plot_w}" y2="{y80}" stroke="#777" stroke-dasharray="5,4"/>')
    for v in [0,.2,.4,.6,.8,1]:
        y=top+(1-v)*plot_h; parts.append(f'<text class="small" x="{left-12}" y="{y+4:.1f}" text-anchor="end">{v:.1f}</text>')
    for p,x in xpos.items(): parts.append(f'<text class="small" x="{x:.1f}" y="{top+plot_h+24}" text-anchor="middle">{p:g}</text>')
    parts.append(f'<text class="lab" x="{left+plot_w/2}" y="{height-55}" text-anchor="middle">Shared-variant prior p12</text>')
    for i,c in enumerate(candidates):
        lx=left+20+i*300; parts.append(f'<circle cx="{lx}" cy="{height-25}" r="5" fill="{colors[c]}"/><text class="small" x="{lx+12}" y="{height-21}">{html.escape(c)}</text>')
    parts.append('</svg>'); path.write_text("\n".join(parts)+"\n")


def main():
    all_rows=list(csv.DictReader(open(ROOT/"coloc_abf_all_priors.tsv"),delimiter="\t"))
    primary=[r for r in all_rows if float(r["p12"])==1e-5]
    windows={r["window_id"]:r for r in csv.DictReader(open(ROOT/"coloc_windows.tsv"),delimiter="\t")}
    fields=["trait","window_id","chrom","start_b37","end_b37","nsnps","duplicate_rows_excluded","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4","H4_over_H3H4","classification"]
    with open(ROOT/"coloc_primary_summary.tsv","w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for r in primary:
            z=windows[r["window_id"]]; w.writerow({**{k:r.get(k) for k in fields},"chrom":z["chrom"],"start_b37":z["start_b37"],"end_b37":z["end_b37"],"classification":classify(r)})
    counts=defaultdict(Counter)
    for r in primary: counts[r["trait"]][classify(r)]+=1
    with open(ROOT/"coloc_classification_counts.tsv","w",newline="") as h:
        w=csv.writer(h,delimiter="\t",lineterminator="\n");w.writerow(["trait","classification","n_windows"])
        for t in counts:
            for c,n in counts[t].items():w.writerow([t,c,n])
    candidates=[f'{r["trait"]}:{r["window_id"]}' for r in primary if classify(r).startswith("PRIMARY")]
    render_overview(primary,ROOT/"coloc_pph4_overview.svg"); render_sensitivity(all_rows,candidates,ROOT/"coloc_prior_sensitivity.svg")


if __name__=="__main__": main()
