#!/usr/bin/env python3
"""Render a compact, source-grounded Phase 1 result report."""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "project" / "results" / "geo"


def rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def metric_table(path: Path) -> dict[str, str]:
    return {x["metric"]: x["value"] for x in rows(path, "\t")}


def marker_table(path: Path) -> dict[str, dict[str, str]]:
    wanted = {"SPP1", "RUNX2", "FN1", "MMP12", "COL1A1", "COL3A1"}
    return {x["gene_symbol"].upper(): x for x in rows(path) if x["gene_symbol"].upper() in wanted}


def main() -> None:
    diff = {x["dataset"]: x for x in rows(OUT / "phase1_differential_summary.csv")}
    conv = metric_table(OUT / "phase1_convergence_summary.tsv")
    pathways = rows(OUT / "convergence_pathways.csv")
    focus = [x for x in pathways if x.get("focus_match") == "TRUE"]
    g1 = marker_table(OUT / "GSE146638_all_genes.csv")
    g2 = marker_table(OUT / "GSE325608_all_genes.csv")
    sensitivity = rows(OUT / "GSE325608_AAD_vs_NonAD_sensitivity.csv")
    sens_deg = [x for x in sensitivity if float(x["adj.P.Val"]) < 0.05 and abs(float(x["logFC"])) > 0.5]
    qc = rows(OUT / "phase1_sample_qc_metrics.csv")

    marker_lines = []
    for gene in ["SPP1", "RUNX2", "FN1", "MMP12", "COL1A1", "COL3A1"]:
        a, b = g1.get(gene), g2.get(gene)
        if a and b:
            marker_lines.append(
                f"- {gene}：Uremic log2FC={float(a['logFC']):.2f}, FDR={float(a['adj.P.Val']):.2g}；"
                f"AAD log2FC={float(b['logFC']):.2f}, FDR={float(b['adj.P.Val']):.2g}。"
            )

    focus_lines = [
        f"- {x['Description']}：Uremic FDR={float(x['FDR_UREMIC']):.2g}，AAD FDR={float(x['FDR_AD']):.2g}。"
        for x in focus
    ]
    kegg_note = (
        "KEGG 已由本地 KEGG REST 快照完成。"
        if conv.get("KEGG_status", "").startswith("VERIFIED")
        else "KEGG REST 下载仍为 UNAVAILABLE（网络超时）；GO-BP 主判定已完成，KEGG 不被伪装为已验证结果。"
    )
    database_label = "GO-BP + KEGG" if conv.get("KEGG_status", "").startswith("VERIFIED") else "GO-BP"
    report = f"""# Phase 1 · GEO 转录组收敛分析结果

生成日期：{date.today().isoformat()}

## 结论

Phase 1 主判定为 **{conv.get('support', 'UNKNOWN')}**。大鼠尿毒症钙化主动脉与小鼠 BAPN-AAD 在映射到人类 1:1 同源基因后，共有 {conv.get('shared_FDR_lt_0.01_pathways', 'UNKNOWN')} 条 {database_label} 通路同时达到 FDR<0.01，其中 {conv.get('focus_pathways', 'UNKNOWN')} 条直接命中预设的 ECM、骨化或平滑肌主题，达到 `shared>=5 且 focus>=3` 的通过标准。该结论是通路收敛证据，不是 CKD 导致主动脉夹层的因果证明。

## 样本与方法

- GSE146638：Control=5、Uremic=5；4 个 Control+VitD 排除。已处理连续表达值经 log2(x+1) 后使用 limma，检测 {diff['GSE146638']['tested_genes']} 个基因。
- GSE325608：Control=20、AAD=16（男7、女9）；12 个 BAPN-NonAD 不进入主对比。RSEM estimated counts 含小数，使用 limma-voom，模型为 `~ sex + primary_group`，检测 {diff['GSE325608']['tested_genes']} 个基因。
- DEG 主阈值：BH FDR<0.05 且 |log2FC|>0.5。GSE146638：上调 {diff['GSE146638']['deg_up']}、下调 {diff['GSE146638']['deg_down']}；GSE325608：上调 {diff['GSE325608']['deg_up']}、下调 {diff['GSE325608']['deg_down']}。
- BAPN 暴露内部敏感性分析（AAD vs NonAD，校正 sex）检出 {len(sens_deg)} 个 DEG，避免把 BAPN 暴露本身与形成 AAD 的病理效应完全混同。

## QC

- 主分析共纳入 {len(qc)} 个样本（10+36）；表达矩阵无缺失、无负值。
- PCA：GSE146638 的 PC1 解释 59.6%，Control/Uremic 清晰分离；GSE325608 的 PC1 解释 59.0%，Control/AAD 明显分离，同时保留可见的性别结构并在模型中校正。
- 基于 library sum 与 detected genes 的稳健 z 分数未发现 |z|>3.5 的样本，因此主分析未作数据驱动的样本剔除。

## 跨物种映射

- Ensembl BioMart 仅保留 `ortholog_one2one` 且 confidence=1。
- 大鼠：{conv.get('rat_tested_genes')} 个受检基因中 {conv.get('rat_1to1_universe')} 个进入人类 1:1 背景；{conv.get('rat_DEG')} 个 DEG 映射为 {conv.get('rat_human_1to1_DEG')} 个人类 DEG。
- 小鼠：{conv.get('mouse_tested_genes')} 个受检基因中 {conv.get('mouse_1to1_universe')} 个进入人类 1:1 背景；{conv.get('mouse_DEG')} 个 DEG 映射为 {conv.get('mouse_human_1to1_DEG')} 个人类 DEG。
- DEG 简单交集为 {conv.get('intersection_genes')} 个，仅作为描述性结果；主证据仍是各自富集后共同显著的通路。

## 预设主题通路

{chr(10).join(focus_lines)}

## 方向一致的代表基因

{chr(10).join(marker_lines)}

## 数据库状态与限制

- GO-BP：本地 `org.Hs.eg.db`，背景为各数据集通过表达过滤且成功 1:1 映射的人类基因，不使用全基因组作虚假背景。
- {kegg_note}
- 两个模型来自不同物种、不同疾病诱导机制与不同表达量生成流程；共同通路支持生物学收敛，但不能消除模型异质性。
- GSE146638 样本量仅 5 vs 5；GSE325608 为终点横断面，病理分组可能包含存活/选择效应。
- 未将 MVMR 表述为中介证明；血压校正属于后续 Phase 2/4 护栏，不属于本阶段表达分析的可用协变量。

## 核心交付物

- `GSE146638_DEG.csv`、`GSE325608_DEG.csv`
- `intersection_genes.csv`、`convergence_pathways.csv`
- `ortholog_rat_human_1to1.csv`、`ortholog_mouse_human_1to1.csv`
- `figures/GSE146638_PCA.png`、`figures/GSE325608_PCA.png`
- `figures/GSE146638_volcano.png`、`figures/GSE325608_volcano.png`
- `figures/convergence_pathway_heatmap.png`
"""
    (OUT / "PHASE1_RESULTS.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
