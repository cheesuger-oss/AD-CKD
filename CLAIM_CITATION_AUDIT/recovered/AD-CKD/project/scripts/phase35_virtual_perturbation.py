#!/usr/bin/env python3
"""Deterministic descriptive perturbation of the frozen undirected STRING network."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, deque
from pathlib import Path

from phase3_network_topology import brandes_betweenness


AD_SET = {"FBN1","ACTA2","COL3A1","COL5A1","EFEMP2","EMILIN1","ENPP1","LOX","MFAP5","MYLK","PLOD3","PRKG1","SMAD2","SMAD3","TGFB2","TGFB3","TGFBR1","TGFBR2","THBS2","THSD4"}
URE_SET = {"ELN","COL1A1","COL3A1","MMP1","MMP2","MMP9","TIMP1","BGN","MGP","RUNX2","SPP1","BMP2","IBSP","NOX4","NFKB1","NFE2L2","ACE","AGTR1","AGT","TGFB1","SMAD3","LOX","TGFBR1","TGFBR2","ENG","BMP4"}
INTERSECTION = sorted(AD_SET & URE_SET)


def components(graph: dict[str, set[str]]) -> list[set[str]]:
    unseen, result = set(graph), []
    while unseen:
        start = min(unseen); seen = {start}; queue = deque([start])
        while queue:
            for neighbor in graph[queue.popleft()]:
                if neighbor not in seen:
                    seen.add(neighbor); queue.append(neighbor)
        unseen -= seen; result.append(seen)
    return result


def distances(graph: dict[str, set[str]], source: str) -> dict[str, int]:
    distance = {source: 0}; queue = deque([source])
    while queue:
        vertex = queue.popleft()
        for neighbor in graph[vertex]:
            if neighbor not in distance:
                distance[neighbor] = distance[vertex] + 1; queue.append(neighbor)
    return distance


def remove_nodes(graph: dict[str, set[str]], removed: set[str]) -> dict[str, set[str]]:
    return {node: neighbors - removed for node, neighbors in graph.items() if node not in removed}


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    source = root / "project/results/network/union_ppi.sif"
    out = root / "project/results/virtual_perturbation"; out.mkdir(parents=True, exist_ok=True)
    graph: dict[str, set[str]] = {}
    with source.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            a, b = row["source"], row["target"]
            graph.setdefault(a, set()).add(b); graph.setdefault(b, set()).add(a)
    missing_ad, missing_ure = sorted(AD_SET - graph.keys()), sorted(URE_SET - graph.keys())
    available_ad = sorted((AD_SET - URE_SET) & graph.keys())
    available_ure = sorted((URE_SET - AD_SET) & graph.keys())
    available_intersection = sorted(set(INTERSECTION) & graph.keys())

    between = brandes_betweenness(graph)
    bet_rows = []
    for gene, value in between.items():
        module = "intersection" if gene in INTERSECTION else ("AD_only" if gene in AD_SET else ("URE_only" if gene in URE_SET else "other"))
        bet_rows.append({"gene": gene, "betweenness_normalized": f"{value:.10f}", "degree": len(graph[gene]), "module": module})
    bet_rows.sort(key=lambda x: (-float(x["betweenness_normalized"]), -int(x["degree"]), str(x["gene"])))
    write_rows(out / "betweenness.csv", bet_rows)

    base_components = components(graph); base_largest = max(map(len, base_components))
    ko_rows = [{"condition": "baseline", "removed_genes": "", "total_nodes": len(graph), "components": len(base_components), "largest_component": base_largest, "largest_fraction": f"{base_largest/len(graph):.8f}", "delta_components": 0, "delta_largest_fraction": "0.00000000"}]
    conditions = [(gene, {gene}) for gene in available_intersection] + [("all_intersection", set(available_intersection))]
    for name, removed in conditions:
        altered = remove_nodes(graph, removed); comps = components(altered); largest = max(map(len, comps))
        fraction = largest / len(altered)
        ko_rows.append({"condition": name, "removed_genes": "|".join(sorted(removed)), "total_nodes": len(altered), "components": len(comps), "largest_component": largest, "largest_fraction": f"{fraction:.8f}", "delta_components": len(comps)-len(base_components), "delta_largest_fraction": f"{fraction-base_largest/len(graph):.8f}"})
    write_rows(out / "ko_fragmentation.csv", ko_rows)

    source_distances = {gene: distances(graph, gene) for gene in available_ure + available_intersection}
    pair_rows, coverage = [], Counter()
    for ure in available_ure:
        for ad in available_ad:
            distance = source_distances[ure].get(ad)
            if distance is None:
                continue
            covered = []
            for gene in available_intersection:
                if gene in source_distances[ure] and ad in source_distances[gene] and source_distances[ure][gene] + source_distances[gene][ad] == distance:
                    covered.append(gene); coverage[gene] += 1
            pair_rows.append({"ure_gene": ure, "ad_gene": ad, "distance": distance, "intersection_on_at_least_one_shortest_path": "|".join(covered)})
    write_rows(out / "ure_ad_pair_distances.csv", pair_rows)
    coverage_rows = [{"gene": gene, "connected_pair_count": len(pair_rows), "pairs_covered": coverage[gene], "percentage": f"{100*coverage[gene]/len(pair_rows):.6f}"} for gene in available_intersection]
    write_rows(out / "path_coverage.csv", coverage_rows)

    summary = {"schema": "ad-ckd.phase35.virtual-perturbation.v1", "graph_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "graph_nodes": len(graph), "graph_edges": sum(map(len, graph.values()))//2, "ad_only_available": len(available_ad), "ure_only_available": len(available_ure), "intersection_available": available_intersection, "missing_ad_seeds": missing_ad, "missing_ure_seeds": missing_ure, "baseline_components": len(base_components), "baseline_largest": base_largest, "ko_all_components": ko_rows[-1]["components"], "ko_all_largest": ko_rows[-1]["largest_component"], "connected_ure_ad_pairs": len(pair_rows), "possible_ure_ad_pairs": len(available_ure)*len(available_ad), "ure_ad_mean_distance": sum(int(row["distance"]) for row in pair_rows)/len(pair_rows), "ure_ad_distance_counts": dict(sorted(Counter(int(row["distance"]) for row in pair_rows).items())), "top_betweenness_gene": bet_rows[0]["gene"], "interpretation": "descriptive robustness of a frozen undirected functional association network; not a biological knockout, causal model, or treatment comparison"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
