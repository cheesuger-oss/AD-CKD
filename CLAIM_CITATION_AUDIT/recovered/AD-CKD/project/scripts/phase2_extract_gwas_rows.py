#!/usr/bin/env python3
"""Stream a gzipped GWAS table and retain rows matching requested rsIDs."""

import argparse
import csv
import gzip
import json
from pathlib import Path


def extract_rows(source: Path, wanted_file: Path, output: Path, snp_column: str) -> dict:
    wanted = {line.strip() for line in wanted_file.read_text().splitlines() if line.strip()}
    rows_scanned = 0
    rows_matched = 0

    with gzip.open(source, "rt") as src:
        header = src.readline().rstrip("\n")
        columns = header.lstrip("#").split("\t")
        if snp_column not in columns:
            raise ValueError(f"SNP column {snp_column!r} is absent from {source}")
        snp_index = columns.index(snp_column)

        with output.open("w") as dst:
            dst.write(header + "\n")
            for line in src:
                rows_scanned += 1
                fields = line.rstrip("\n").split("\t")
                if snp_index >= len(fields):
                    continue
                row_ids = set(fields[snp_index].split(","))
                if row_ids & wanted:
                    dst.write(line)
                    rows_matched += 1

    return {
        "source": str(source),
        "wanted": len(wanted),
        "rows_scanned": rows_scanned,
        "rows_matched": rows_matched,
        "output": str(output),
    }


def extract_rows_by_positions(
    source: Path, wanted_file: Path, output: Path, chromosome_column: str, position_column: str
) -> dict:
    with wanted_file.open(newline="") as handle:
        wanted = {
            (row["chr"], row["pos"])
            for row in csv.DictReader(handle, delimiter="\t")
        }
    rows_scanned = rows_matched = 0
    with gzip.open(source, "rt") as src:
        header = src.readline().rstrip("\n")
        columns = header.lstrip("#").split("\t")
        for column in (chromosome_column, position_column):
            if column not in columns:
                raise ValueError(f"Position column {column!r} is absent from {source}")
        chromosome_index = columns.index(chromosome_column)
        position_index = columns.index(position_column)
        with output.open("w") as dst:
            dst.write(header + "\n")
            for line in src:
                rows_scanned += 1
                fields = line.rstrip("\n").split("\t")
                if max(chromosome_index, position_index) >= len(fields):
                    continue
                if (fields[chromosome_index], fields[position_index]) in wanted:
                    dst.write(line)
                    rows_matched += 1
    return {
        "source": str(source),
        "wanted": len(wanted),
        "rows_scanned": rows_scanned,
        "rows_matched": rows_matched,
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("wanted_file", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("snp_column")
    parser.add_argument("--position-columns", nargs=2, metavar=("CHR", "POS"))
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    if args.position_columns:
        summary = extract_rows_by_positions(
            args.source, args.wanted_file, args.output, *args.position_columns
        )
    else:
        summary = extract_rows(args.source, args.wanted_file, args.output, args.snp_column)
    if args.summary:
        args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
