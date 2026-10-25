# Corrected revision code and Source Data — author-approved public archive

CBAC-D-26-05362; public-file scope approved on 25 October 2026.
This revision archive is separate from the original registered Zenodo version
1.0.0 (doi:10.5281/zenodo.21937216); that record does not establish coverage of
this revision. The approved repository is https://github.com/cheesuger-oss/AD-CKD.
Repository publication and download verification are recorded separately from
file preparation; the manuscript will identify the verified revision commit.

## Which results are current
CLAIM_CITATION_AUDIT/identifier_correction_20261024 contains the current corrected
rat results and plot inputs. Earlier recovered/candidate inputs are historical
dependencies, not substitutes for the corrected values. The corrected workbook
is retained at UPLOAD_READY_R1_FINAL6/SuppTables_R1.xlsx because the original S3
wrapper uses that path; the directory name does not confer release approval.

25 date-converted rat identifiers were resolved using original annotation and
rn5 coordinates. Only affected rat workflows were recomputed with the original
samples, models and thresholds. Tested rat genes: 9814; mapped background: 7292;
DEGs: 1734; shared animal ORA pathways: 154. Shared genes remain 549, comprising
478 concordant and 71 discordant. Six frozen focus pathways were not reselected.
Rat KEGG vascular smooth muscle contraction: NES -1.56, BH q=0.0555; negative
enrichment but not BH-significant. Mouse, human, genetic and classifier analyses
were not refitted. The GSTM2 table identity correction is distinct from this
rat-identifier recalculation.

## Display regeneration from supplied results
Run from the extracted root, using R 4.5 on macOS (Quartz, ggplot2, patchwork,
scales) and Python with matplotlib/openpyxl. Session records are included.
The following commands render plots only; they do not refit statistical models:

    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/render_corrected_core_figures.R
    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/render_corrected_f4_f5.R
    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/render_corrected_s8.R
    python CLAIM_CITATION_AUDIT/identifier_correction_20261024/render_corrected_s3.py
    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/render_ga_journal_size.R

Statistical plot outputs go to identifier_correction_20261024/figures. The current
graphical abstract goes to identifier_correction_20261024/journal_figure_candidate/figures.
The historical render_corrected_ga.R is retained for provenance, not the current
graphical-abstract export command. The S3 generator includes positional labels
so that dataset identification does not depend on colour alone. Never replace submission
files without reviewing exports. Full raw-data statistical reproduction requires
public source inputs, rn5/NCBI reference files, targeted orthology responses and
the frozen annotation environment; this is not a self-contained full-rerun claim.
Included recalculation scripts document the workflow and require the documented
external dependencies. Identifier verification takes the original workbook path
from ADCKD_GSE146638_RAW_XLSX. Historical table generators require an explicit
ADCKD_REBUILD_HISTORICAL=1 opt-in and accept ADCKD_PROJECT_ROOT; they produce
historical tables, not the current reviewed tables. Do not run them as display steps.
Historical network directory names do not denote biological intervention.

MANIFEST.csv records payload sizes and SHA-256 hashes; SHA256SUMS.txt also covers
this README and manifest. Hashes are integrity checks, not watermarks. Third-party
data retain their original terms. No private correspondence or credentials are
included. Isolated export execution and complete package review remain separate
gates; a valid archive checksum alone is not approval.

## Author-approved public scope

The author approved the exact 500-file public scope. See LICENSING.md for the
author-approved MIT and CC BY 4.0 scope and third-party exclusions. The
current supplied workbooks take precedence over historical table generators.
No statistical outputs were changed during path sanitization.

## Restricted-data exclusions

See PUBLIC_MATRIX_EXCLUSIONS.json and THIRD_PARTY_NOTICES.md. Three expression
matrix files are intentionally excluded; original inputs must be reacquired
from GEO for affected analysis scripts. This is not a complete raw-data rerun
archive. The full local audit package is separate and is not included in the
public upload. Approval of the public-file scope is not a claim that all original
raw-data analyses have been independently reproduced.
