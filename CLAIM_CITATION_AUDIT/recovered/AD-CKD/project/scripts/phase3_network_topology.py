#!/usr/bin/env python3
"""Compute deterministic descriptive topology for the frozen undirected STRING network."""

from __future__ import annotations

import csv
import json
from collections import deque
from pathlib import Path


def brandes_betweenness(adjacency: dict[str, set[str]]) -> dict[str, float]:
    """Normalized undirected Brandes betweenness, matching common graph-library convention."""
    nodes = sorted(adjacency)
    score = dict.fromkeys(nodes, 0.0)
    for source in nodes:
        stack: list[str] = []
        predecessors = {node: [] for node in nodes}
        paths = dict.fromkeys(nodes, 0.0)
        paths[source] = 1.0
        distance = dict.fromkeys(nodes, -1)
        distance[source] = 0
        queue = deque([source])
        while queue:
            vertex = queue.popleft()
            stack.append(vertex)
            for neighbor in sorted(adjacency[vertex]):
                if distance[neighbor] < 0:
                    queue.append(neighbor)
                    distance[neighbor] = distance[vertex] + 1
                if distance[neighbor] == distance[vertex] + 1:
                    paths[neighbor] += paths[vertex]
                    predecessors[neighbor].append(vertex)
        dependency = dict.fromkeys(nodes, 0.0)
        while stack:
            child = stack.pop()
            if paths[child]:
                coefficient = (1.0 + dependency[child]) / paths[child]
                for parent in predecessors[child]:
                    dependency[parent] += paths[parent] * coefficient
            if child != source:
                score[child] += dependency[child]
    n = len(nodes)
    scale = 1.0 / ((n - 1) * (n - 2)) if n > 2 else 0.0
    return {node: value * scale for node, value in score.items()}


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    base = root / "project/results/network"
    adjacency: dict[str, set[str]] = {}
    weighted_degree: dict[str, float] = {}
    with (base / "union_ppi.sif").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            a, b, weight = row["source"], row["target"], float(row["score"])
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
            weighted_degree[a] = weighted_degree.get(a, 0.0) + weight
            weighted_degree[b] = weighted_degree.get(b, 0.0) + weight
    between = brandes_betweenness(adjacency)
    rows = [
        {
            "gene": node, "degree": len(adjacency[node]),
            "weighted_degree": f"{weighted_degree[node]:.6f}",
            "betweenness_normalized": f"{between[node]:.10f}",
            "rank_rule": "degree_desc;betweenness_desc;gene_asc",
        }
        for node in adjacency
    ]
    rows.sort(key=lambda row: (-int(row["degree"]), -float(row["betweenness_normalized"]), str(row["gene"])))
    with (base / "hub_genes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema": "ad-ckd.phase3.topology.v1", "node_count": len(adjacency),
        "edge_count": sum(len(value) for value in adjacency.values()) // 2,
        "ranking_rule": "degree descending, then normalized betweenness descending, then gene ascending",
        "betweenness_definition": "unweighted undirected normalized Brandes",
        "weighted_degree_definition": "sum of STRING functional association scores",
        "interpretation": "descriptive network topology only; hub rank is not causal evidence",
    }
    (base / "network_topology_metadata.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
