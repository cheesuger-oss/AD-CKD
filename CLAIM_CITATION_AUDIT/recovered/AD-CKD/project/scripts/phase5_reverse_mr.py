#!/usr/bin/env python3
"""Prespecified reverse MR from two independent AAD instruments."""

import csv, gzip, math
from pathlib import Path


TRAITS=["egfrcrea","egfrcys","cystatin_c","bun","ckd","uacr"]
LEADS={"rs17352842","rs1232604"}


def normal_p(z): return math.erfc(abs(z)/math.sqrt(2))


def ivw(bx,by,sey):
    w=[1/(s*s) for s in sey]; den=sum(wi*x*x for wi,x in zip(w,bx)); num=sum(wi*x*y for wi,x,y in zip(w,bx,by))
    b=num/den; se=math.sqrt(1/den); q=sum(wi*(y-b*x)**2 for wi,x,y in zip(w,bx,by))
    # df=1 here, so chi-square upper tail equals erfc(sqrt(Q/2)).
    qp=math.erfc(math.sqrt(q/2))
    return b,se,q,qp


def bh_adjust(p):
    n=len(p); order=sorted(range(n),key=p.__getitem__); out=[0.0]*n; running=1.0
    for rank0,i in reversed(list(enumerate(order))):
        running=min(running,p[i]*n/(rank0+1));out[i]=running
    return out


def main():
    root=Path("project/results/shared_genetics/coloc"); results=[]; instruments=[]
    for trait in TRAITS:
        found={}
        with gzip.open(root/f"aad__{trait}_harmonized.tsv.gz","rt") as h:
            for r in csv.DictReader(h,delimiter="\t"):
                if r["rsid"] in LEADS: found[r["rsid"]]=r
        for r in found.values():
            instruments.append({"trait":trait,"rsid":r["rsid"],"beta_aad":r["aad_beta"],"se_aad":math.sqrt(float(r["aad_varbeta"])),"beta_outcome":r["renal_beta"],"se_outcome":math.sqrt(float(r["renal_varbeta"]))})
        if len(found)<2:
            results.append({"trait":trait,"n_instruments":len(found),"beta":None,"se":None,"ci_low":None,"ci_high":None,"p":None,"q_bh":None,"Q":None,"Q_p":None,"status":"UNDERIDENTIFIED"});continue
        rows=[found[x] for x in sorted(LEADS)]; bx=[float(r["aad_beta"]) for r in rows];by=[float(r["renal_beta"]) for r in rows];sey=[math.sqrt(float(r["renal_varbeta"])) for r in rows]
        b,se,q,qp=ivw(bx,by,sey);results.append({"trait":trait,"n_instruments":2,"beta":b,"se":se,"ci_low":b-1.96*se,"ci_high":b+1.96*se,"p":normal_p(b/se),"q_bh":None,"Q":q,"Q_p":qp,"status":"ESTIMABLE"})
    estimable=[r for r in results if r["status"]=="ESTIMABLE"]
    for r,q in zip(estimable,bh_adjust([r["p"] for r in estimable])):r["q_bh"]=q
    for name,rows in [("reverse_mr_instruments.tsv",instruments),("reverse_mr_results.tsv",results)]:
        with open(root/name,"w",newline="") as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    est=[r for r in results if r["status"]=="ESTIMABLE"]; width,height,left,top=820,410,160,45; xlo=min(r["ci_low"] for r in est);xhi=max(r["ci_high"] for r in est);pad=(xhi-xlo)*.08;xlo-=pad;xhi+=pad
    sx=lambda x:left+(x-xlo)/(xhi-xlo)*600
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><rect width="100%" height="100%" fill="white"/>','<style>text{font-family:Arial,sans-serif;fill:#222}.lab{font-size:15px}.small{font-size:12px}</style>',f'<line x1="{sx(0):.1f}" y1="25" x2="{sx(0):.1f}" y2="340" stroke="#777" stroke-dasharray="4,4"/>']
    for i,r in enumerate(est):
        y=top+i*58;parts += [f'<text class="lab" x="{left-15}" y="{y+5}" text-anchor="end">{r["trait"]}</text>',f'<line x1="{sx(r["ci_low"]):.1f}" y1="{y}" x2="{sx(r["ci_high"]):.1f}" y2="{y}" stroke="#0072B2" stroke-width="3"/>',f'<circle cx="{sx(r["beta"]):.1f}" cy="{y}" r="5" fill="#0072B2"/>',f'<text class="small" x="{sx(r["ci_high"])+8:.1f}" y="{y+4}">q={r["q_bh"]:.3f}</text>']
    parts += [f'<text class="lab" x="{left+300}" y="385" text-anchor="middle">Renal-trait effect per unit genetically predicted AAD log-odds (95% CI)</text>','</svg>']
    (root/"reverse_mr_forest.svg").write_text("\n".join(parts)+"\n")
    for r in results:print(r)


if __name__=="__main__":main()
