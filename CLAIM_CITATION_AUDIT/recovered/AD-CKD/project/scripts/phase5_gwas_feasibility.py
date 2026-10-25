#!/usr/bin/env python3
"""Stream full GWAS files and freeze feasibility facts for shared-genetics analyses."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import TextIO


P_ALIASES = ("pval", "p_value", "p-value", "p.value.gc", "p-value.gc", "p")
RSID_ALIASES = ("rsids", "rsid", "hm_rsid", "snp", "markername")


def _split(line: str, delimiter: str) -> list[str]:
    return line.rstrip("\r\n").split("\t") if delimiter == "tab" else line.split()


def _find_column(header: list[str], aliases: tuple[str, ...]) -> int | None:
    normalized = {name.lstrip("#").casefold(): i for i, name in enumerate(header)}
    for alias in aliases:
        if alias.casefold() in normalized:
            return normalized[alias.casefold()]
    return None


def audit_stream(handle: TextIO, delimiter: str) -> dict:
    header_line = handle.readline()
    if not header_line:
        raise ValueError("empty GWAS stream")
    header = _split(header_line, delimiter)
    p_idx = _find_column(header, P_ALIASES)
    rsid_idx = _find_column(header, RSID_ALIASES)
    result = {
        "columns": header,
        "p_column": header[p_idx] if p_idx is not None else None,
        "rsid_column": header[rsid_idx] if rsid_idx is not None else None,
        "variant_rows": 0,
        "valid_p_rows": 0,
        "gws_rows": 0,
        "rsid_rows": 0,
        "min_p": None,
    }
    for raw in handle:
        if not raw.strip():
            continue
        fields = _split(raw, delimiter)
        result["variant_rows"] += 1
        if rsid_idx is not None and rsid_idx < len(fields):
            value = fields[rsid_idx].strip()
            if value and value not in {".", "NA", "NaN", "nan"}:
                result["rsid_rows"] += 1
        if p_idx is None or p_idx >= len(fields):
            continue
        try:
            p = float(fields[p_idx])
        except (TypeError, ValueError):
            continue
        if not 0 <= p <= 1:
            continue
        result["valid_p_rows"] += 1
        if p < 5e-8:
            result["gws_rows"] += 1
        if result["min_p"] is None or p < result["min_p"]:
            result["min_p"] = p
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_file(path: Path, delimiter: str, include_hash: bool = False) -> dict:
    opener = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        result = audit_stream(handle, delimiter)
    result.update({"path": str(path), "bytes": path.stat().st_size})
    if include_hash:
        result["sha256"] = sha256(path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--delimiter", choices=("tab", "whitespace"), required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hash", action="store_true")
    args = parser.parse_args()
    result = audit_file(Path(args.file), args.delimiter, args.hash)
    result["label"] = args.label
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("label", "variant_rows", "gws_rows", "min_p")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
