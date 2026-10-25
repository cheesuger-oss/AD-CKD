"""Query only corrected identifiers; retain original mapping for every other gene."""
from pathlib import Path
import json,urllib.request,urllib.parse,hashlib
R=Path(__file__).resolve().parent
v=json.loads((R/'identity_verification.json').read_text());assert v['status']=='PASS'
symbols=sorted({r['corrected_symbol'] for r in v['records']})
xml='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="rnorvegicus_gene_ensembl" interface="default"><Filter name="external_gene_name" value="'+','.join(symbols)+'"/>'+''.join('<Attribute name="'+a+'"/>' for a in ['external_gene_name','ensembl_gene_id','hsapiens_homolog_associated_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_orthology_type','hsapiens_homolog_orthology_confidence'])+'</Dataset></Query>'
(R/'reference/targeted_biomart_query.xml').write_text(xml)
url='https://www.ensembl.org/biomart/martservice?'+urllib.parse.urlencode({'query':xml})
with urllib.request.urlopen(url,timeout=60) as f:data=f.read()
(R/'reference/targeted_biomart_response.tsv').write_bytes(data)
print(data.decode()[:9000])
