# Independent human expression validation specification

Frozen: 2026-08-09, before inspecting pathway-validation results in the selected datasets.

## Objective

Test whether the six Phase 1 focus pathways reproduce in independent human aortic-dissection tissue and human uremia-associated arterial tissue. This is molecular-convergence validation, not evidence that CKD causes incident aortic dissection.

## Dataset-selection rules

Datasets were screened using species, disease definition, tissue, platform, group size, sample-level expression availability and obvious group-batch confounding before testing the six pathways. Phase 1 genes, OMIM/uremia seeds, STRING hubs, DGIdb targets and genetic results were not used to choose datasets.

### Primary AAD dataset

- GSE153434: human ascending aorta, Stanford type A aortic dissection versus non-aortic-disease CABG controls, bulk RNA-seq, 10 versus 10.
- Primary contrast: TAAD minus normal.
- Raw integer counts supplied by GEO will be analysed with edgeR/limma-voom after expression filtering and library normalization.
- Known limitation: surgical controls differ clinically from emergency TAAD patients; metadata do not support adjustment for age, sex or other covariates.

### Primary renal-vascular dataset

- GSE38750 (SubSeries of GSE38752): human arterial tissue, uremic renal-transplant recipients versus healthy kidney donors, Affymetrix HG-U133A 2.0, 15 versus 19.
- Primary contrast: uremia minus control.
- The GEO normalized series matrix will be analysed with limma after probe-to-Entrez mapping; when several probes map to one gene, the probe with the largest across-sample variance is retained without reference to case-control effect.
- Critical design limitation: case tissue is iliac artery and control tissue is renal artery, so disease is perfectly confounded with arterial bed. Results are supportive/hypothesis-generating and cannot isolate a pure uremic effect.

## Exclusions and secondary datasets

- GSE160611, GSE212388 and GSE266437 are flow/shear experiments in human aortic endothelial cells; `UF` means unidirectional flow, not uremic filtrate. They are not CKD validation datasets.
- GSE189795, GSE213740, GSE222318, GSE226492 and GSE254132 are human AAD single-cell datasets. They may later localize signals by cell type but will not be pooled with bulk samples or pseudoreplicated at cell level.
- GSE147026, GSE190635, GSE52093 and GSE98770 are eligible independent human AAD bulk datasets for sensitivity/replication after the primary analysis.
- GSE146638 and GSE325608 are the Phase 1 animal discovery datasets and cannot count as human replication.

## Frozen pathway family

The primary family contains exactly six GO Biological Process terms:

1. GO:0030198 extracellular matrix organization
2. GO:0001503 ossification
3. GO:0006939 smooth muscle contraction
4. GO:0048659 smooth muscle cell proliferation
5. GO:0048660 regulation of smooth muscle cell proliferation
6. GO:0048661 positive regulation of smooth muscle cell proliferation

Membership will come from the same local human GO annotation framework used in Phase 1. It will not be redefined using genes significant in the validation datasets.

## Analysis and quality gates

- Confirm sample order, group sizes, numeric range, missingness, duplicated identifiers, library/distribution behaviour and gene mapping before modelling.
- Produce PCA and sample-correlation heatmaps for every primary dataset. No sample will be removed solely to improve separation; any exclusion requires a prespecified technical failure or an explicit sensitivity analysis.
- Differential-expression tables will report effect, standard error/test statistic, raw P and BH FDR for all tested genes. Volcano and top-variable/differential heatmap figures are descriptive QC/supporting figures, not the primary validation test.
- Primary pathway test: preranked GSEA using the full gene-level moderated t-statistic, without DEG prefiltering. Within each dataset, BH correction is applied across the six frozen pathways.
- A pathway is supported within a dataset when BH q<0.05. Direction and normalized enrichment score are reported; nested smooth-muscle-proliferation terms are not described as independent biological replications.
- The human-convergence hypothesis is supported only if at least one technically suitable human dataset on each disease side shows corrected enrichment within the frozen family. Because GSE38750 has tissue-bed confounding, even a positive result remains supportive rather than definitive.

## Interpretation prohibitions

- No threshold, dataset or contrast will be changed in response to non-significant results.
- Human expression enrichment does not establish direction, mediation or incident disease risk.
- Phase 2 MVMR is not mediation proof, and network hubs cannot be used to infer causality or upgrade expression evidence.

