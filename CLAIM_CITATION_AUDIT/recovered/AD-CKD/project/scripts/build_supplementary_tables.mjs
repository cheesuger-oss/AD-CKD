import fs from "node:fs/promises";
import path from "node:path";
import os from "node:os";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

// Historical generator: its old result counts must not overwrite revision tables.
if (process.env.ADCKD_REBUILD_HISTORICAL !== "1") throw new Error("Historical table generator only. Current reviewed tables are supplied separately. Set ADCKD_REBUILD_HISTORICAL=1 only to reproduce the old table with its original inputs.");
const root = path.resolve(process.env.ADCKD_PROJECT_ROOT || path.join(path.dirname(fileURLToPath(import.meta.url)), "../.."));
const output = path.join(root, "project/manuscript/outputs/ad_ckd_submission/Supplementary_Tables_S1-S4.xlsx");

const specs = [
  ["S1", "S1_rat_DEG", "project/results/geo/GSE146638_DEG.csv"],
  ["S1", "S1_rat_orth_DEG", "project/results/geo/GSE146638_human_1to1_DEG.csv"],
  ["S1", "S1_mouse_DEG", "project/results/geo/GSE325608_DEG.csv"],
  ["S1", "S1_mouse_orth_DEG", "project/results/geo/GSE325608_human_1to1_DEG.csv"],
  ["S2", "S2_GO_AAD", "project/results/geo/GOBP_AD.csv"],
  ["S2", "S2_GO_uremic", "project/results/geo/GOBP_UREMIC.csv"],
  ["S2", "S2_KEGG_AAD", "project/results/geo/KEGG_AD.csv"],
  ["S2", "S2_KEGG_uremic", "project/results/geo/KEGG_UREMIC.csv"],
  ["S2", "S2_GSEA_all", "project/results/geo/GSEA_all_results.csv"],
  ["S2", "S2_convergence", "project/results/geo/convergence_pathways.csv"],
  ["S3", "S3_topology", "project/results/network/hub_genes.csv"],
  ["S3", "S3_enrichment", "project/results/network/enrichment_table.csv"],
  ["S3", "S3_betweenness", "project/results/virtual_perturbation/betweenness.csv"],
  ["S3", "S3_node_removal", "project/results/virtual_perturbation/ko_fragmentation.csv"],
  ["S3", "S3_seed_distances", "project/results/virtual_perturbation/ure_ad_pair_distances.csv"],
  ["S3", "S3_path_coverage", "project/results/virtual_perturbation/path_coverage.csv"],
  ["S4", "S4_immune", "project/results/phase4/analysis/immune_program_contrasts.tsv"],
  ["S4", "S4_random_controls", "project/results/phase4/analysis/random_gene_control_auc.tsv"],
  ["S4", "S4_external_audit", "project/results/phase4/analysis/external_validation_results.tsv"],
  ["S4", "S4_stability", "project/results/phase4/analysis/feature_stability.tsv"],
  ["S4", "S4_coefficients", "project/results/phase4/analysis/locked_signature_coefficients.tsv"],
];

function parseDelimited(text, delimiter) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === delimiter) { row.push(field); field = ""; }
    else if (c === "\n") { row.push(field.replace(/\r$/, "")); rows.push(row); row = []; field = ""; }
    else field += c;
  }
  if (field.length || row.length) { row.push(field.replace(/\r$/, "")); rows.push(row); }
  return rows;
}

function typedRows(rows) {
  if (!rows.length) return rows;
  const headers = rows[0];
  const keepText = headers.map(h => /(^|_)(id|gene|dataset|method|module|condition|status|description|outcome|exposure|rank_rule)($|_)/i.test(h));
  return rows.map((row, ri) => row.map((v, ci) => {
    if (ri === 0 || keepText[ci] || v === "") return v;
    if (/^(TRUE|FALSE)$/i.test(v)) return v.toUpperCase() === "TRUE";
    if (/^-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(v)) return Number(v);
    return v;
  }));
}

function colName(n) {
  let s = "";
  for (n++; n; n = Math.floor((n - 1) / 26)) s = String.fromCharCode(65 + ((n - 1) % 26)) + s;
  return s;
}

const wb = Workbook.create();
const readme = wb.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:D1").values = [["Supplement", "Worksheet", "Frozen source file", "Interpretation boundary"]];
const indexRows = specs.map(([group, sheet, source]) => [
  `Supplementary Table ${group}`,
  sheet,
  source,
  group === "S3" ? "Descriptive database-network topology; no causal or therapeutic inference" :
  group === "S4" ? "Immune programme scores are not cell proportions; retrospective tissue classification is not future-risk prediction" :
  group === "S2" ? "ORA/GSEA are pathway-level enrichment results, not mediation evidence" :
  "Complete gene-level results under the frozen contrasts and mapping rules"
]);
readme.getRange(`A2:D${indexRows.length + 1}`).values = indexRows;
readme.getRange("A1:D1").format = { fill: "#1F4E79", font: { bold: true, color: "#FFFFFF" }, wrapText: true };
readme.getRange(`A1:D${indexRows.length + 1}`).format.borders = { preset: "inside", style: "thin", color: "#D9E2F3" };
readme.getRange("A:A").format.columnWidth = 22;
readme.getRange("B:B").format.columnWidth = 24;
readme.getRange("C:C").format.columnWidth = 68;
readme.getRange("D:D").format.columnWidth = 74;
readme.getRange(`A2:D${indexRows.length + 1}`).format.wrapText = true;
readme.freezePanes.freezeRows(1);

for (let i = 0; i < specs.length; i++) {
  const [, sheetName, rel] = specs[i];
  const delimiter = rel.endsWith(".tsv") ? "\t" : ",";
  const text = await fs.readFile(path.join(root, rel), "utf8");
  let rows = typedRows(parseDelimited(text.replace(/^\uFEFF/, ""), delimiter));
  if (rel.endsWith("external_validation_results.tsv")) {
    const statusCol = rows[0].indexOf("status");
    rows[0][statusCol] = "source_status";
    rows[0].splice(statusCol + 1, 0, "reporting_status");
    for (let r = 1; r < rows.length; r++) {
      const reporting = rows[r][statusCol] === "NOT_READY" ? "NOT EVALUABLE" : rows[r][statusCol];
      rows[r].splice(statusCol + 1, 0, reporting);
    }
  }
  const width = rows[0].length;
  if (!rows.every(r => r.length === width)) throw new Error(`Non-rectangular source: ${rel}`);
  const sheet = wb.worksheets.add(sheetName);
  sheet.showGridLines = false;
  const range = sheet.getRange(`A1:${colName(width - 1)}${rows.length}`);
  range.values = rows;
  sheet.getRange(`A1:${colName(width - 1)}1`).format = {
    fill: "#1F4E79", font: { bold: true, color: "#FFFFFF" }, wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#17365D" }
  };
  sheet.getRange(`A2:${colName(width - 1)}${Math.min(rows.length, 500)}`).format.borders = {
    preset: "inside", style: "thin", color: "#E7E6E6"
  };
  range.format.autofitColumns();
  for (let c = 0; c < width; c++) {
    const col = sheet.getRange(`${colName(c)}:${colName(c)}`);
    if (col.format.columnWidth > 38) col.format.columnWidth = 38;
  }
  sheet.freezePanes.freezeRows(1);
  sheet.tables.add(`A1:${colName(width - 1)}${rows.length}`, true, `T${String(i + 1).padStart(2, "0")}`);
}

await fs.mkdir(path.dirname(output), { recursive: true });
const out = await SpreadsheetFile.exportXlsx(wb);
await out.save(output);

const inspect = await wb.inspect({ kind: "sheet", include: "id,name", maxChars: 5000 });
console.log(inspect.ndjson);
const errors = await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "formula error scan" });
console.log(errors.ndjson);
for (const name of ["README", ...specs.map(x => x[1])]) {
  const png = await wb.render({ sheetName: name, range: name === "README" ? "A1:D12" : "A1:H12", scale: 1 });
  await fs.writeFile(path.join(os.tmpdir(), `ad_ckd_${name}.png`), new Uint8Array(await png.arrayBuffer()));
}
console.log(`OUTPUT=${output}`);
