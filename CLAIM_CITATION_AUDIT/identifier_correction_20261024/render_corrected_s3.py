"""Correction-only inputs for the original Python S3 plotting workflow."""
from pathlib import Path
import csv,hashlib,json,runpy,os
from matplotlib.colorbar import Colorbar
root=Path(__file__).resolve().parent
audit=root.parent
frozen=audit/'supplementary_source_inputs/AD-CKD/project/results/geo'
paths={
 'S2_GO_uremic':root/'results/corrected_GO_BP_UREMIC.csv',
 'S2_KEGG_uremic':root/'results/corrected_KEGG_UREMIC.csv',
 'S2_GO_AAD':frozen/'GOBP_AD.csv',
 'S2_KEGG_AAD':frozen/'KEGG_AD.csv',
 'old_gsea':frozen/'GSEA_all_results.csv',
 'rat_go_gsea':root/'results/corrected_GSEA_GO_BP.csv',
 'rat_kegg_gsea':root/'results/corrected_GSEA_KEGG.csv'}
numeric={'p.adjust','pvalue','qvalue','FoldEnrichment','Count','NES','setSize','RichFactor','zScore','enrichmentScore','rank'}
def read(path):
 with path.open() as f: data=list(csv.DictReader(f))
 for r in data:
  for k in numeric & r.keys(): r[k]=None if r[k] in ('','NA') else float(r[k])
 return data
def rows(sheet):
 if sheet!='S2_GSEA_all':return read(paths[sheet])
 data=[r for r in read(paths['old_gsea']) if r['dataset']!='Uremic']
 for key,db in [('rat_go_gsea','GO_BP'),('rat_kegg_gsea','KEGG')]:
  for r in read(paths[key]):
   r.update(dataset='Uremic',database=db);data.append(r)
 return data
original=Colorbar._draw_all
def vector_draw(self):
 original(self)
 if self.solids is not None:
  self.solids.set_rasterized(False);self.solids.set_edgecolor('face')
Colorbar._draw_all=vector_draw
source=audit/'assemble_supplementary_figure_s3_AAD.py'
code=source.read_text()
old_label='Blue: uraemic calcification; red: AAD lesions'
assert code.count(old_label)==1
code=code.replace(old_label,'Lower ten: uraemic calcification; upper ten: AAD lesions')
ns={'__file__':str(source),'__name__':'s3_accessible_source'}
exec(compile(code,str(source),'exec'),ns)
out=Path(os.environ.get('CBC_FIGURE_OUTPUT',str(root/'figures')));out.mkdir(exist_ok=True,parents=True)
ns['main'].__globals__.update(ROOT=out,rows=rows,
 OUT_PDF=out/'Supplementary_Figure_S3_animal_enrichment_summary.pdf',
 OUT_PNG=out/'Supplementary_Figure_S3_animal_enrichment_summary.png')
ns['main']()
path=out/'S3_publication_layout_source.json'
payload=json.loads(path.read_text())
payload.pop('source_workbook_sha256',None)
payload.update(source_inputs={k:{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},
 new_analysis=False,upstream_recalculation='Identifier correction; rat-only affected analyses recomputed with original models and thresholds',
 selection_rule_unchanged=True)
path.write_text(json.dumps(payload,indent=2))
