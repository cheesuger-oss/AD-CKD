"""Export S3 panels by executing the original Python drawing code, not PDF cropping."""
from pathlib import Path
import json,subprocess,hashlib
r=Path(__file__).resolve().parent
out=r/'independent_panels'
out.mkdir(exist_ok=True)
ns={'__file__':str(r/'render_corrected_s3.py'),'__name__':'panel_export'}
wrapper=(r/'render_corrected_s3.py').read_text()
# Load the established corrected frozen-result reader and vector colorbar hook.
exec(compile(wrapper.split('source=audit/')[0],str(r/'render_corrected_s3.py'),'exec'),ns)
read_frozen=ns['rows']
source=r.parent/'assemble_supplementary_figure_s3_AAD.py'
code=source.read_text().replace('Blue: uraemic calcification; red: AAD lesions','Lower ten: uraemic calcification; upper ten: AAD lesions')
code=code.replace('fontsize=8.0, fontweight="bold")','fontsize=8.0, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=2))')
code=code.replace('fontsize=7.5, color=GREY)','fontsize=7.5, color=GREY, bbox=dict(facecolor="white", edgecolor="none", pad=2))')
start=code.index('    fig = plt.figure(');end=code.index('    # A: five',start)
code=code[:start]+'''    figures=[]
    axes=[]
    for panel,height in zip('ABCDEF',[155,155,95,95,95,95]):
        f,ax=plt.subplots(figsize=(183/25.4,height/25.4),constrained_layout=True)
        figures.append(f);axes.append(ax)
    fig=figures[0];ax_a=axes[0];ax_b=axes[1];mini=axes[2:]
'''+code[end:]
start=code.index('    fig.suptitle(');end=code.index("    payload=",start)
code=code[:start]+'''    for label,f in zip('ABCDEF',figures):
        f.savefig(ROOT/('FigureS3_'+label+'.pdf'))
        plt.close(f)
'''+code[end:]
exec(compile(code,str(source),'exec'),ns)
ns.update(rows=read_frozen,ROOT=out,BOOK=r/'documents/SuppTables_R1.xlsx',OUT_PDF=out/'FigureS3_A.pdf',OUT_PNG=out/'FigureS3_A.png')
ns['main']()
made=json.loads((out/'S3_publication_layout_source.json').read_text())
frozen=json.loads((r/'journal_figure_candidate/figures/S3_publication_layout_source.json').read_text())
assert all(made[k]==frozen[k] for k in ('ora','go','gsea_selected'))
made.pop('page_heights_mm',None)
made['panel_heights_mm']=dict(zip('ABCDEF',[155,155,95,95,95,95]))
made['export_mode']='independent panels from original drawing code'
(out/'S3_publication_layout_source.json').write_text(json.dumps(made,indent=2)+'\n')
for label in 'ABCDEF':
    pdf=out/f'FigureS3_{label}.pdf'
    subprocess.run(['/opt/homebrew/bin/pdftoppm','-scale-to','1700','-singlefile','-png',str(pdf),str(pdf.with_suffix(''))],check=True)
(out/'s3_panel_export.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'corrected_reader_sha256':hashlib.sha256((r/'render_corrected_s3.py').read_bytes()).hexdigest(),'panels':list('ABCDEF'),'ORA_GO_GSEA_selected_rows_identical':True,'visual_review':'PENDING','new_analysis':False},indent=2)+'\n')
