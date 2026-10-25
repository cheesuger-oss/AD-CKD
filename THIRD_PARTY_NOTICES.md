# Third-party sources and rights — proposed release

Checked 24 October 2026. This notice supplements, and does not expand, LICENSING.md.
This notice is included in the local candidate ZIP. No upload is authorized by its inclusion.
MIT applies only to author-owned code. CC BY 4.0 applies only to contributions the
AD-CKD authors have authority to license. Neither license overrides source-provider rights.

## GEO expression data

Expression-derived matrices and sample-level diagnostic values originate from the
accessions identified in the source tables and manuscript. NCBI itself does not
restrict use or distribution of GEO data, but expressly preserves potential
submitter intellectual-property rights. Public accessibility is not a declaration
that the authors own the underlying measurements. Retain GEO accession identifiers
and the original study citations, including GSE146638, GSE325608, GSE153434,
GSE38750 and GSE38751. Corrections to rat identifiers and derived calculations are
documented separately; the original investigators did not endorse those changes.

Source: https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html
Dataset-specific restrictions must still be checked where an actual expression
matrix is redistributed; this notice does not certify their absence.

## Gene Ontology

Gene Ontology data and data products are provided by the Gene Ontology Consortium
under CC BY 4.0. Preserve attribution to the Consortium, the original source link,
license and disclaimer. The analysis-used annotation/package versions, not the
current website release, govern reproducibility. Do not silently replace them.
Our enrichment statistics are newly calculated results, not a new GO release.

Source and disclaimer: https://www.geneontology.org/page/use-and-license
License: https://creativecommons.org/licenses/by/4.0/

## STRING

STRING data and outputs are CC BY 4.0, with attribution and identification of changes.
The study uses the frozen human STRING v12.0 functional association network;
topology and path metrics were calculated by the AD-CKD authors. Material originally
from other providers retains its own terms. This archive must not suggest STRING
endorsement of the study or of a causal/therapeutic interpretation.

Source: https://en.string-db.org/cgi/access?footer_active_subpage=licensing
License: https://creativecommons.org/licenses/by/4.0/

## KEGG and MSigDB

KEGG is copyright Kanehisa Laboratories and is not blanket-relicensed under this
project's CC BY 4.0 notice. The candidate contains enrichment results (statistics,
pathway identifiers/names and observed contributing genes), rather than a full KEGG
or MSigDB gene-set database. No permission to redistribute a complete KEGG database
is asserted. The published analysis results and original source identifiers remain
attributed; users should obtain reference databases from their original providers
and comply with the corresponding terms.

MSigDB releases have collection-specific exceptions: KEGG MEDICUS and KEGG LEGACY
must not be treated as having an identical license. No MSigDB GMT/XML database is
included in the inspected candidate inventory. Do not add one without a new review.

Sources: https://www.genome.jp/kegg/legal.html
https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp

## Ensembl and orthology

Ensembl-generated data are available without restriction, while incorporated
third-party content may have separate constraints. Preserve original orthology
provenance and the analysis-used versions. The corrected rat-to-human mapping is
an analysis input/derivative, not a claim of ownership over external annotations.

Source: https://www.ensembl.org/info/about/legal/disclaimer.html

## GWAS resources and software dependencies

The inspected candidate contains study-level MR, colocalization and LDSC outputs
and a source registry, not full downloadable GWAS summary-statistic releases or
individual-level biobank data. Obtain original GWAS inputs through their providers;
access agreements and acknowledgements still apply. FinnGen results should credit
the FinnGen participants and investigators and the source publication.

Source: https://www.finngen.fi/en/access_results

R/Python packages and external software are dependencies, not relicensed by MIT
in this project. Preserve any notices in reused code. Historical table generators
are marked as historical and do not define the current revision results.

## Remaining release conditions

This is a source-rights assessment, not an unrestricted legal warranty. Complete
dataset-specific checks and file-level attribution before final publication.
The final ZIP must include the approved notice and regenerated checksums, and the
author must confirm the final inventory before upload.
