#!/usr/bin/env python3
import csv, gzip, json, math, sys
from pathlib import Path

src, out, state = map(Path, sys.argv[1:4])
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 500000
done = json.loads(state.read_text()) if state.exists() else {"rows": 0, "offset": 0, "hits": 0, "complete": False}
keep = ["chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR"]
with gzip.open(src, "rb") as f:
    header = f.readline().decode().rstrip().split("\t")
    idx = {x: header.index(x) for x in keep + ["low_confidence_EUR"]}
    if done["offset"]: f.seek(done["offset"])
    mode = "a" if out.exists() else "w"
    with out.open(mode, newline="") as g:
        w = csv.writer(g, delimiter="\t")
        if mode == "w": w.writerow(keep)
        n = hits = 0
        while n < limit:
            raw = f.readline()
            if not raw: break
            row = raw.decode().rstrip().split("\t"); n += 1
            try:
                b, se, lp = float(row[idx["beta_EUR"]]), float(row[idx["se_EUR"]]), float(row[idx["neglog10_pval_EUR"]])
                if lp > 7.30103 and se > 0 and math.isfinite(b) and math.isfinite(se) and row[idx["low_confidence_EUR"]].lower() == "false":
                    w.writerow([row[idx[x]] for x in keep]); hits += 1
            except (ValueError, IndexError): pass
        done.update(rows=done["rows"]+n, offset=f.tell(), hits=done["hits"]+hits, complete=n < limit)
state.write_text(json.dumps(done)+"\n")
print(json.dumps(done), flush=True)
