"""Retain frozen mappings and append only coordinate-verified corrected identifiers."""
from pathlib import Path
import csv,gzip,json
R=Path(__file__).resolve().parent
def read(p,sep=','):
    with p.open() as f:return list(csv.DictReader(f,delimiter=sep))
identity=json.loads((R/'identity_verification.json').read_text());assert identity['status']=='PASS'
with gzip.open(R/'reference/Rattus_norvegicus.gene_info.gz','rt') as f:genes={x['GeneID']:x for x in csv.DictReader(f,delimiter='\t')}
bio=read(R/'reference/targeted_biomart_response.tsv','\t')
old=read(R.parent/'supplementary_source_inputs/AD-CKD/project/results/geo/ortholog_rat_human_1to1.csv')
new=list(old);audit=[]
for record in identity['records']:
    gene=genes[record['ncbi_gene_ids'][0]]
    ids=[x.split(':',1)[1] for x in gene['dbXrefs'].split('|') if x.startswith('Ensembl:')]
    hits=[x for x in bio if x['Gene stable ID'] in ids]
    accepted=[x for x in hits if x['Human homology type']=='ortholog_one2one' and x['Human orthology confidence [0 low, 1 high]']=='1' and x['Human gene name']]
    for x in accepted:
        new.append({'source_symbol':record['corrected_symbol'],'human_symbol':x['Human gene name'],'orthology_type':x['Human homology type'],'orthology_confidence':'1'})
    audit.append({'corrected_symbol':record['corrected_symbol'],'ncbi_gene_id':gene['GeneID'],'stable_ids':ids,'biomart_rows':hits,'accepted_human_symbols':[x['Human gene name'] for x in accepted],
                  'join':'Stable Ensembl ID anchored by NCBI cross-reference; not a symbol-only join'})
assert len({(x['source_symbol'],x['human_symbol']) for x in new})==len(new)
with (R/'ortholog_rat_human_1to1_corrected.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(old[0]));w.writeheader();w.writerows(new)
oldDEG=read(R/'results/baseline_DEG.csv');newDEG=read(R/'results/corrected_DEG.csv')
def mapped(rows,maps):
    d={x['gene_symbol']:x for x in rows};return {x['human_symbol']:('up' if float(d[x['source_symbol']]['logFC'])>0 else 'down') for x in maps if x['source_symbol'] in d}
oldCandidates=mapped(oldDEG,old);newCandidates=mapped(newDEG,new)
out={'added_mapping_rows':len(new)-len(old),'corrected_gene_mapping_audit':audit,
     'renal_human_DEG_count_before':len(oldCandidates),'renal_human_DEG_count_after':len(newCandidates),
     'renal_candidate_gene_and_direction_sets_identical':oldCandidates==newCandidates,
     'unchanged_mapping_rows_preserved':new[:len(old)]==old,
     'mapping_scope':'Targeted current Ensembl lookup for the 25 repaired identities only; frozen mappings elsewhere unchanged. Apply original high-confidence one-to-one criterion. Do not force Septin4 into human mapping when no qualifying result is supplied.'}
(R/'mapping_dependency_audit.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='corrected_gene_mapping_audit'})
