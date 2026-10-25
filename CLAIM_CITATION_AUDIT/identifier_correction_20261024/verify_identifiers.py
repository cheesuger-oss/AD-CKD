"""Verify source annotation against rn5 RefSeq loci and NCBI rat identifier aliases."""
from pathlib import Path
import csv,gzip,datetime,json,hashlib,os
from openpyxl import load_workbook
R=Path(__file__).resolve().parent
raw_path=os.environ.get('ADCKD_GSE146638_RAW_XLSX')
if not raw_path:
    raise SystemExit('Set ADCKD_GSE146638_RAW_XLSX to the original GEO processed workbook; this verification also requires the documented reference files.')
raw=Path(raw_path)
ref=[]
with gzip.open(R/'reference/rn5_refGene.txt.gz','rt') as f:
    for row in csv.reader(f,delimiter='\t'):ref.append({'chrom':row[2],'start0':int(row[4]),'end0':int(row[5]),'symbol':row[12],'transcript':row[1]})
aliases={}
with gzip.open(R/'reference/Rattus_norvegicus.gene_info.gz','rt') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        for key in [row['Symbol'],*row['Synonyms'].split('|')]:
            if key!='-':aliases.setdefault(key,[]).append(row)
w=load_workbook(raw,read_only=True,data_only=True)
records=[]
for i,row in enumerate(w.active.values,1):
    if not isinstance(row[0],datetime.datetime):continue
    chrom,interval=row[2].split(':');start,end=map(int,interval.split('-'))
    hits=[r for r in ref if r['chrom']==chrom and min(end,r['end0'])>max(start,r['start0'])]
    # Source locus convention is not assumed; candidates must remain identical under a 1-bp start shift.
    assert {h['symbol'] for h in hits}=={r['symbol'] for r in ref if r['chrom']==chrom and min(end,r['end0'])>max(start-1,r['start0'])}
    name=row[1]
    expected=('Sept'+name.split()[-1]) if name.startswith('septin ') else ('Sep15' if name=='selenoprotein 15' else ('Marc'+name.split()[-1] if name.startswith('mitochondrial amidoxime') else 'March'+name.split('finger')[1].split(',')[0].replace('(C3HC4)','').strip()))
    matched=[]
    for h in hits:
        for gene in aliases.get(h['symbol'],[]):
            if expected in [gene['Symbol'],*gene['Synonyms'].split('|')] or name.casefold()==gene['description'].casefold():matched.append((h,gene))
    ids={g['GeneID'] for _,g in matched}
    record={'source_row':i,'date_identifier':str(row[0]),'source_annotation':name,'source_locus_rn5':row[2],
            'annotation_legacy_symbol':expected,'overlapping_refseq_symbols':sorted({h['symbol'] for h in hits}),
            'ncbi_gene_ids':sorted(ids),'verified':len(ids)==1,
            'corrected_symbol':matched[0][1]['Symbol'] if len(ids)==1 else None,
            'refseq_transcripts':sorted({h['transcript'] for h,g in matched})}
    records.append(record);print(i,name,record['overlapping_refseq_symbols'],record['corrected_symbol'],flush=True)
report={'assembly':'rn5 / Rnor_5.0 per GEO processing description',
 'coordinate_policy':'RefSeq genePred is 0-based half-open; source locus convention not explicitly documented; overlapping candidate sets verified invariant under 1-base start shift. Identity additionally requires annotation-derived alias and unique NCBI Gene ID.',
 'reference_urls':['https://hgdownload.soe.ucsc.edu/goldenPath/rn5/database/refGene.txt.gz','https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4403082'],
 'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw,*sorted((R/'reference').glob('*.gz'))] if not p.name.startswith('._')},
 'records':records,'status':'PASS' if len(records)==25 and all(x['verified'] for x in records) else 'REVIEW_REQUIRED'}
(R/'identity_verification.json').write_text(json.dumps(report,indent=2));print(report['status'])
if report['status']=='PASS':
    mapping={x['source_row']:x for x in records}
    original=R.parent/'identifier_source_audit/GSE146638_expression.tsv.gz'
    with gzip.open(original,'rt') as f:previous=list(csv.reader(f,delimiter='\t'))
    source_rows=list(w.active.values)
    assert len(previous)==len(source_rows)
    changes=0
    output=[]
    for i,(src,old) in enumerate(zip(source_rows,previous),1):
        if i==1:
            assert old[1:]==list(src[12:]);output.append(old);continue
        assert old[0]==str(src[0])
        assert all(float(a)==float(b) for a,b in zip(old[1:],src[12:])),i
        new=old.copy()
        if i in mapping:new[0]=mapping[i]['corrected_symbol'];changes+=1
        output.append(new)
    assert changes==25
    with gzip.open(R/'GSE146638_expression_corrected.tsv.gz','wt',newline='') as f:csv.writer(f,delimiter='\t').writerows(output)
    (R/'expression_preservation.json').write_text(json.dumps({'changed_identifiers':changes,'expression_values_unchanged':True,'rows':len(output)-1,'samples':len(output[0])-1,'original_gzip_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'corrected_gzip_sha256':hashlib.sha256((R/'GSE146638_expression_corrected.tsv.gz').read_bytes()).hexdigest()},indent=2))
