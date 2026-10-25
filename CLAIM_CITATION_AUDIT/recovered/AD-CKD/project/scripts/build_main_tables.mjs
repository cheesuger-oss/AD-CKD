import fs from "node:fs/promises";
import path from "node:path";
import os from "node:os";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

// Historical generator: its old result counts must not overwrite revision tables.
if (process.env.ADCKD_REBUILD_HISTORICAL !== "1") throw new Error("Historical table generator only. Current reviewed tables are supplied separately. Set ADCKD_REBUILD_HISTORICAL=1 only to reproduce the old table with its original inputs.");
const root = path.resolve(process.env.ADCKD_PROJECT_ROOT || path.join(path.dirname(fileURLToPath(import.meta.url)), "../.."));
const output = path.join(root, "project/manuscript/outputs/ad_ckd_submission/Main_Tables_1-2.xlsx");

function parseDelimited(text, delimiter) {
  const rows=[]; let row=[], field="", quoted=false;
  for (let i=0;i<text.length;i++) {
    const c=text[i];
    if (quoted) { if(c==='"'&&text[i+1]==='"'){field+='"';i++;} else if(c==='"')quoted=false; else field+=c; }
    else if(c==='"') quoted=true;
    else if(c===delimiter){row.push(field);field="";}
    else if(c==='\n'){row.push(field.replace(/\r$/,"") );rows.push(row);row=[];field="";}
    else field+=c;
  }
  if(field.length||row.length){row.push(field.replace(/\r$/,"") );rows.push(row);}
  return rows;
}
function objects(rows){const h=rows[0];return rows.slice(1).map(r=>Object.fromEntries(h.map((x,i)=>[x,r[i]??""])));}
async function load(rel,delim){return objects(parseDelimited((await fs.readFile(path.join(root,rel),"utf8")).replace(/^\uFEFF/,""),delim));}
function n3(x){return Number(x).toFixed(3);}

const ratMeta=await load("project/results/geo/GSE146638_sample_metadata.tsv","\t");
const mouseMeta=await load("project/results/geo/GSE325608_sample_metadata.tsv","\t");
const ratDeg=await load("project/results/geo/GSE146638_DEG.csv",",");
const mouseDeg=await load("project/results/geo/GSE325608_DEG.csv",",");
const ratOrth=await load("project/results/geo/GSE146638_human_1to1_DEG.csv",",");
const mouseOrth=await load("project/results/geo/GSE325608_human_1to1_DEG.csv",",");
const renal=await load("project/results/mr/renal_family/renal_family_mr_results.tsv","\t");
const mainMr=await load("project/results/mr/main_mr_results.csv",",");

const count=(xs,key,val)=>xs.filter(x=>x[key]===val).length;
const direction=(xs,up)=>xs.filter(x=>Number(x["adj.P.Val"])<0.05 && (up?Number(x.logFC)>0.5:Number(x.logFC)<-0.5)).length;

const table1=[
 ["Feature","Uremic aortic model","Aortic dissection model"],
 ["Accession","GSE146638","GSE325608"],
 ["Species",ratMeta[0].species,mouseMeta[0].species],
 ["Primary groups",`Control, n=${count(ratMeta,"primary_group","Control")}; uremic, n=${count(ratMeta,"primary_group","Uremic")}`,`Control, n=${count(mouseMeta,"primary_group","Control")}; AAD, n=${count(mouseMeta,"primary_group","AAD")}`],
 ["Excluded from primary contrast",`Control + vitamin D, n=${ratMeta.filter(x=>x.include_primary==="no").length}`,`BAPN-exposed non-AAD, n=${mouseMeta.filter(x=>x.include_primary==="no").length}`],
 ["Expression input","Processed continuous expression","RSEM estimated counts containing non-integer values"],
 ["Primary model","limma after log2(x+1)","limma-voom, adjusted for sex"],
 ["Tested genes",9813,14326],
 ["DEG threshold","BH-adjusted P<0.05 and |log2FC|>0.5","BH-adjusted P<0.05 and |log2FC|>0.5"],
 ["Upregulated / downregulated DEG",`${direction(ratDeg,true)} / ${direction(ratDeg,false)}`,`${direction(mouseDeg,true)} / ${direction(mouseDeg,false)}`],
 ["Human 1:1 orthologue background",7275,12391],
 ["Human 1:1 orthologous DEG",ratOrth.length,mouseOrth.length],
];

const names={egfrcrea:"eGFRcrea",serum_cystatin_c:"Serum cystatin C",egfrcys:"eGFRcys",bun_ea:"Blood urea nitrogen",ckd_ea:"Chronic kidney disease",microalbuminuria:"Microalbuminuria",rapid3:"Rapid3",ckdi25:"CKDi25",uacr_ea:"Urinary albumin-to-creatinine ratio"};
const primary=[...mainMr.filter(x=>x.method==="Inverse variance weighted").map(x=>({...x,exposure_id:"egfrcrea",q_primary_family:"0.685111596064231"})),...renal.filter(x=>x.primary_ivw==="TRUE")];
const order=["egfrcrea","serum_cystatin_c","egfrcys","bun_ea","ckd_ea","microalbuminuria","rapid3","ckdi25","uacr_ea"];
const table2=[["Renal phenotype","Harmonized variants","Primary method","Odds ratio (95% CI)","P value","BH q value"]];
for(const id of order){const x=primary.find(r=>r.exposure_id===id);table2.push([names[id],Number(x.nsnp),"IVW",`${n3(x.OR)} (${n3(x.OR_lci95)}–${n3(x.OR_uci95)})`,Number(x.pval),Number(x.q_primary_family)]);}

const wb=Workbook.create();
for(const [name,title,data,note] of [
 ["Table 1","Table 1 | Transcriptomic datasets and frozen primary analyses",table1,"The two human orthologous DEG sets shared 549 genes. AAD, aortic dissection; BAPN, β-aminopropionitrile; BH, Benjamini–Hochberg; DEG, differentially expressed gene."],
 ["Table 2","Table 2 | Primary renal phenotype Mendelian-randomization estimates for aortic dissection",table2,"BH correction was applied jointly across the nine prespecified primary renal-phenotype estimates. None passed the family false-discovery threshold; this is not evidence of no effect. CI, confidence interval; IVW, inverse-variance weighted."]]){
  const s=wb.worksheets.add(name);s.showGridLines=false;
  const endCol=name==="Table 1"?"C":"F"; const start=3,end=2+data.length;
  s.getRange(`A1:${endCol}1`).merge();s.getRange("A1").values=[[title]];
  s.getRange(`A1:${endCol}1`).format={fill:"#17365D",font:{bold:true,color:"#FFFFFF",size:14},rowHeight:28};
  s.getRange(`A${start}:${endCol}${end}`).values=data;
  s.getRange(`A${start}:${endCol}${start}`).format={fill:"#1F4E79",font:{bold:true,color:"#FFFFFF"},wrapText:true,borders:{preset:"outside",style:"thin",color:"#17365D"}};
  s.getRange(`A${start+1}:${endCol}${end}`).format.borders={preset:"inside",style:"thin",color:"#D9E2F3"};
  s.getRange(`A${start+1}:${endCol}${end}`).format.wrapText=true;
  s.tables.add(`A${start}:${endCol}${end}`,true,name==="Table 1"?"MainTable1":"MainTable2");
  const noteRow=end+2;s.getRange(`A${noteRow}:${endCol}${noteRow+1}`).merge();s.getRange(`A${noteRow}`).values=[[note]];s.getRange(`A${noteRow}`).format={wrapText:true,font:{italic:true,color:"#404040"}};
  s.getRange(`A1:${endCol}${noteRow+1}`).format.autofitColumns();
  for(let c=0;c<(name==="Table 1"?3:6);c++){const col=String.fromCharCode(65+c);if(s.getRange(`${col}:${col}`).format.columnWidth>42)s.getRange(`${col}:${col}`).format.columnWidth=42;}
  if(name==="Table 1") { s.getRange("B:C").format.columnWidth=48; s.getRange(`A3:C${end}`).format.autofitRows(); }
  if(name==="Table 2"){s.getRange(`E4:F${end}`).format.numberFormat="0.000";}
  s.freezePanes.freezeRows(3);
}
const p=wb.worksheets.add("Provenance");p.showGridLines=false;
const provenance=[
 ["Table","Frozen source","Use"],
 ["Table 1","project/results/geo/GSE146638_sample_metadata.tsv; project/results/geo/GSE325608_sample_metadata.tsv","Sample counts and exclusions"],
 ["Table 1","project/results/geo/GSE146638_DEG.csv; project/results/geo/GSE325608_DEG.csv","DEG counts and directions"],
 ["Table 1","project/results/geo/GSE146638_human_1to1_DEG.csv; project/results/geo/GSE325608_human_1to1_DEG.csv","Mapped DEG counts"],
 ["Table 1","project/manuscript/METHODS_MAIN_APPROVED_SECTIONS.md","Frozen tested-gene/background definitions and analysis models"],
 ["Table 2","project/results/mr/main_mr_results.csv; project/results/mr/renal_family/renal_family_mr_results.tsv","Nine primary IVW estimates and family BH q values"],
];
p.getRange(`A1:C${provenance.length}`).values=provenance;p.getRange("A1:C1").format={fill:"#1F4E79",font:{bold:true,color:"#FFFFFF"}};p.getRange(`A1:C${provenance.length}`).format.wrapText=true;p.getRange("A:A").format.columnWidth=14;p.getRange("B:B").format.columnWidth=80;p.getRange("C:C").format.columnWidth=42;p.freezePanes.freezeRows(1);

await fs.mkdir(path.dirname(output),{recursive:true});const out=await SpreadsheetFile.exportXlsx(wb);await out.save(output);
console.log((await wb.inspect({kind:"sheet",include:"id,name",maxChars:2000})).ndjson);
console.log((await wb.inspect({kind:"match",searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",options:{useRegex:true,maxResults:100},summary:"formula error scan"})).ndjson);
for(const [name,range] of [["Table 1","A1:C18"],["Table 2","A1:F16"],["Provenance","A1:C7"]]){const png=await wb.render({sheetName:name,range,scale:2});await fs.writeFile(path.join(os.tmpdir(), `ad_ckd_main_${name.replace(/ /g,"_")}.png`),new Uint8Array(await png.arrayBuffer()));}
console.log(`OUTPUT=${output}`);
