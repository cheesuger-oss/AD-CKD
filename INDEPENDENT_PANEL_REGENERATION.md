# Independent panel regeneration

Work on a disposable extracted copy of this archive: the original display
scripts also write composite outputs. No statistical model is refitted.

The tested graphics backend is macOS R with Quartz and ggplot2/patchwork;
S3 retains its original Python/matplotlib backend. Poppler must be installed
at /opt/homebrew/bin for the existing preview commands. Other platforms are
not yet tested; this is not a platform-independent environment lockfile.

From the extracted Revision_Code_SourceData_R1 directory:

    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/export_independent_core_panels.R CLAIM_CITATION_AUDIT/identifier_correction_20261024 core
    Rscript CLAIM_CITATION_AUDIT/identifier_correction_20261024/export_independent_core_panels.R CLAIM_CITATION_AUDIT/identifier_correction_20261024 remaining
    python3 CLAIM_CITATION_AUDIT/identifier_correction_20261024/export_independent_s3.py

Install the R/Python packages recorded in the supplied session/dependency
information. Python S3 requires matplotlib, numpy and openpyxl. No package
installation occurs automatically. Inspect executable scripts before running.
These commands use supplied derived plot inputs; they do not require the
three excluded full expression matrices. Full analysis reproduction is a
different workflow and requires the original public inputs.

GSE146638 originals: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146638
Use the study sample workbook and its original sample identities, not a new
cohort. The correction scripts document identifier repair and duplicate
collapse. Do not substitute a different database release or gene universe.
See PUBLIC_MATRIX_EXCLUSIONS.json for exact omitted files and checksums.

Status: candidate pending isolated regeneration test and final release review.
