"""Stage only affected panels; retain unrelated rows byte-for-value without regeneration."""
from pathlib import Path
import csv,json,math,hashlib
R=Path(__file__).resolve().parent;S=R/'results';F=R.parent/'aad_revision_candidate/figures';O=R/'corrected_figure_source_data';O.mkdir(exist_ok=True)
def read(p,sep='\t'):
    with p.open() as f:return list(csv.DictReader(f,delimiter=sep))
changes=[]
def write(name,rows):
    with (O/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    changes.append({'file':name,'rows_before':len(read(F/name)),'rows_after':len(rows),'before_sha256':hashlib.sha256((F/name).read_bytes()).hexdigest(),'after_sha256':hashlib.sha256((O/name).read_bytes()).hexdigest()})
rat=read(S/'corrected_all_genes.csv',',')
name='Figure1_animal_transcriptomic_discovery_source_data_B_renal_volcano.tsv'
write(name,[{k:r[k] for k in ['gene_symbol','logFC','adj.P.Val']}|{'status':('Up' if float(r['logFC'])>0 else 'Down') if float(r['adj.P.Val'])<.05 and abs(float(r['logFC']))>.5 else 'Not DEG'} for r in rat])
ora={r['ID']:r for r in read(S/'corrected_GO_BP_UREMIC.csv',',')}
gsea={r['ID']:r for r in read(S/'corrected_GSEA_GO_BP.csv',',')}
for suffix in ['A_ORA','B_animal_GSEA']:
    name=f'Figure2_pathway_convergence_human_support_source_data_{suffix}.tsv';rows=read(F/name)
    for r in rows:
        if r['model']!='Renal dysfunction':continue
        if suffix=='A_ORA':r['q']=ora[r['pathway_id']]['p.adjust']
        else:
            for k in ['NES','pvalue','p.adjust']:r[k]=gsea[r['pathway_id']][k]
            r['passed']='TRUE' if float(r['p.adjust'])<.05 else 'FALSE'
    write(name,rows)
for prefix,suffix,key,numerics in [
 ('Supplementary_Figure_S1_animal_QC','A_PCA','PCA',['PC1','PC2']),
 ('Supplementary_Figure_S1_animal_QC','C_correlation','correlation',['rho']),
 ('Supplementary_Figure_S2_animal_expression_detail','A_heatmap','heatmap',['row_z'])]:
    name=f'{prefix}_source_data_{suffix}.tsv';old=read(F/name);base=read(S/f'baseline_{key}.tsv');new=read(S/f'corrected_{key}.tsv')
    assert len(old)==len(base)
    assert all(abs(float(a[k])-float(b[k]))<1e-9 for a,b in zip(old,base) for k in numerics),(name,'baseline failure')
    # Keep the established column order. All sample rows come from frozen metadata.
    write(name,[{k:r[k] for k in old[0]} for r in new])
name='Figure4_immune_classifier_robustness_source_data_A_immune.tsv';rows=read(F/name)
new={r['program_id']:r for r in read(S/'corrected_immune_contrasts.tsv')}
for r in rows:
    if r['dataset']!='GSE146638':continue
    for k,v in new[r['program_id']].items():r[k]=v
    effect=float(r['standardized_effect']);n1=int(r['case_n']);n0=int(r['control_n'])
    se=math.sqrt((n1+n0)/(n1*n0)+effect**2/(2*(n1+n0-2)))
    r.update(effect_se=str(se),ci_low=str(effect-1.96*se),ci_high=str(effect+1.96*se))
write(name,rows)
name='Supplementary_Figure_S8_source_data_A_immune_scores.tsv';rows=read(F/name)
new={(r['sample'],r['program_id']):r for r in read(S/'corrected_immune_scores_display.tsv')}
for r in rows:
    if r['dataset']=='GSE146638':r.update(new[r['sample'],r['program_id']])
write(name,rows)
report={'scope':'Corrected panel source data staged separately; not yet promoted into submission package','files':changes,'baseline_QC_reproduction':'PASS','human_mouse_genetic_rows':'not recomputed'}
(R/'figure_source_data_change_manifest.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
