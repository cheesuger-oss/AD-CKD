#!/usr/bin/env python3
"""Render the LDSC QC forest plot as a dependency-free SVG."""

import csv
import html
import sys


def main(input_tsv, output_svg):
    with open(input_tsv, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    width, height = 920, 120 + 62 * len(rows)
    left, right, top = 230, 40, 55
    plot_w = width - left - right
    xmax = max(float(r["ci_high"]) for r in rows) * 1.08
    x = lambda value: left + float(value) / xmax * plot_w
    colors = {"PASS": "#0072B2", "CAUTION": "#E69F00", "BLOCK": "#D55E00"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.lab{font-size:16px}.small{font-size:13px}</style>']
    axis_y = top + 62 * len(rows) + 10
    parts.append(f'<line x1="{left}" y1="{top-20}" x2="{left}" y2="{axis_y}" stroke="#777" stroke-dasharray="3,4"/>')
    for i, row in enumerate(rows):
        y = top + i * 62
        color = colors[row["qc_status"]]
        parts.append(f'<text class="lab" x="{left-15}" y="{y+5}" text-anchor="end">{html.escape(row["trait"])}</text>')
        parts.append(f'<line x1="{x(row["ci_low"]):.1f}" y1="{y}" x2="{x(row["ci_high"]):.1f}" y2="{y}" stroke="{color}" stroke-width="4"/>')
        parts.append(f'<circle cx="{x(row["h2_observed"]):.1f}" cy="{y}" r="6" fill="{color}"/>')
        parts.append(f'<text class="small" x="{x(row["ci_high"])+9:.1f}" y="{y+5}">{row["h2_observed"]} ({row["qc_status"]})</text>')
    parts.append(f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>')
    for tick in [0, 0.1, 0.2, 0.3, 0.4]:
        if tick <= xmax:
            tx = x(tick)
            parts.extend([f'<line x1="{tx:.1f}" y1="{axis_y}" x2="{tx:.1f}" y2="{axis_y+6}" stroke="#222"/>',
                          f'<text class="small" x="{tx:.1f}" y="{axis_y+23}" text-anchor="middle">{tick:.1f}</text>'])
    parts.append(f'<text class="lab" x="{left+plot_w/2:.1f}" y="{height-16}" text-anchor="middle">Observed-scale SNP heritability (95% CI)</text>')
    parts.append('</svg>')
    with open(output_svg, "w") as handle:
        handle.write("\n".join(parts) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: plot_ldsc_h2_qc.py input.tsv output.svg")
    main(sys.argv[1], sys.argv[2])
