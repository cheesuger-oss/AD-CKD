#!/usr/bin/env python3
import csv, gzip, json, sys
from pathlib import Path

src, out, state = map(Path, sys.argv[1:4])
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 750000
done = json.loads(state.read_text()) if state.exists() else {"rows": 0, "offset": 0, "header": None, "complete": False}
if done["complete"]:
    print(json.dumps(done)); raise SystemExit(0)
with gzip.open(src, "rb") as f:
    header = f.readline().decode().rstrip("\r\n").split("\t")
    if done.get("offset", 0): f.seek(done["offset"])
    elif done["rows"]:
        for _ in range(done["rows"]): f.readline()
    mode = "a" if out.exists() else "w"
    with out.open(mode, newline="") as g:
        writer = csv.writer(g, delimiter="\t")
        if mode == "w": writer.writerow(header)
        n = 0
        while True:
            if n >= limit: break
            raw = f.readline()
            if not raw: break
            row = raw.decode().rstrip("\r\n").split("\t")
            n += 1
            try:
                if float(row[header.index("P.value.GC")]) < 5e-8 and row[header.index("RSID")].startswith("rs"):
                    writer.writerow(row)
            except (ValueError, IndexError):
                pass
        done["rows"] += n
        done["offset"] = f.tell()
        done["header"] = header
        done["complete"] = n < limit
state.write_text(json.dumps(done, ensure_ascii=False) + "\n")
print(json.dumps(done), flush=True)
