#!/usr/bin/env python3
import gzip, json, sys
from pathlib import Path

src, ids_file, out, state = map(Path, sys.argv[1:5])
limit = int(sys.argv[5]) if len(sys.argv) > 5 else 500000
wanted = {x.strip() for x in ids_file.read_text().splitlines() if x.strip()}
done = json.loads(state.read_text()) if state.exists() else {"rows": 0, "offset": 0, "hits": 0, "complete": False}
with gzip.open(src, "rb") as f:
    header = f.readline().decode()
    if done["offset"]: f.seek(done["offset"])
    mode = "a" if out.exists() else "w"
    with out.open(mode) as g:
        if mode == "w": g.write(header)
        n = hits = 0
        while n < limit:
            raw = f.readline()
            if not raw: break
            line = raw.decode(); n += 1
            fields = line.split("\t", 6)
            if len(fields) > 4 and fields[4] in wanted:
                g.write(line); hits += 1
        done["rows"] += n; done["hits"] += hits; done["offset"] = f.tell(); done["complete"] = n < limit
state.write_text(json.dumps(done) + "\n")
print(json.dumps(done), flush=True)
