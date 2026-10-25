#!/usr/bin/env Rscript

# Programmatic, non-quantitative synthesis of the study's evidential narrowing.
# Box dimensions are deliberately constant and do not encode effect magnitude.

suppressPackageStartupMessages(library(grid))

out_dir <- "project/results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(out_dir, "Figure5_evidence_contraction.png")
pdf_file <- file.path(out_dir, "Figure5_evidence_contraction.pdf")
tsv_file <- file.path(out_dir, "Figure5_evidence_contraction_source_data.tsv")

stages <- data.frame(
  order = 1:5,
  stage = c(
    "Initial observation",
    "Human tissue test",
    "Genetic test",
    "Specificity and transfer tests",
    "Retained claim"
  ),
  evidence = c(
    "Animal aortas: 549 shared DEGs; 155 shared pathways; matrix and smooth-muscle programmes",
    "Two of six frozen programmes recurred in both datasets: extracellular matrix organization and ossification",
    "No primary MR estimate passed BH correction; LDSC non-estimable; colocalization prior-sensitive; reverse MR limited",
    "No immune programme passed correction in either human dataset; random-gene controls performed similarly; external AUCs varied",
    "A bounded arterial-remodelling tissue-state overlap, strongest for matrix-related programmes"
  ),
  interpretation = c(
    "Broad molecular convergence in experimental models",
    "Narrower, non-directional human support; one dataset confounds disease with arterial site",
    "Does not establish a renal-to-dissection relationship; smaller effects remain possible",
    "Does not establish immune replication, renal-origin specificity or future-risk prediction",
    "Not a causal mechanism, therapeutic target or generalizable renal-derived signature"
  ),
  status = c("Observed", "Narrowed", "Constrained", "Constrained", "Supported boundary"),
  stringsAsFactors = FALSE
)
write.table(stages, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cols <- c(
  "Initial observation" = "#0072B2",
  "Human tissue test" = "#D55E00",
  "Genetic test" = "#CC79A7",
  "Specificity and transfer tests" = "#6B7280",
  "Retained claim" = "#009E73"
)

wrap <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 8, ncol = 1,
    heights = unit(c(13, 4, 19, 19, 19, 19, 24, 7), "mm")
  )))

  grid.text(
    "Evidence progressively narrows the initial molecular convergence",
    vp = viewport(layout.pos.row = 1), x = unit(0.02, "npc"), just = "left",
    gp = gpar(fontfamily = "sans", fontsize = 12.5, fontface = "bold", col = "#111827")
  )

  for (i in 1:4) {
    row <- i + 2
    col <- unname(cols[stages$stage[i]])
    pushViewport(viewport(layout.pos.row = row, x = 0.5, width = 0.96))
    grid.roundrect(gp = gpar(fill = "#FFFFFF", col = "#D1D5DB", lwd = 0.9), r = unit(2, "mm"))
    grid.rect(x = 0.012, width = 0.024, just = "left", gp = gpar(fill = col, col = NA))
    grid.text(stages$stage[i], x = 0.045, y = 0.70, just = "left",
              gp = gpar(fontfamily = "sans", fontsize = 8.3, fontface = "bold", col = col))
    grid.text(wrap(stages$evidence[i], 91), x = 0.045, y = 0.38, just = "left",
              gp = gpar(fontfamily = "sans", fontsize = 7.2, col = "#111827", lineheight = 1.08))
    grid.text(wrap(stages$interpretation[i], 60), x = 0.985, y = 0.50, just = "right",
              gp = gpar(fontfamily = "sans", fontsize = 6.9, col = "#4B5563", lineheight = 1.08))
    popViewport()
    if (i < 4) {
      grid.lines(x = unit(c(0.5, 0.5), "npc"),
                 y = unit(c((6 - i) / 8 + 0.005, (6 - i) / 8 - 0.005), "npc"),
                 arrow = arrow(type = "closed", length = unit(1.6, "mm")),
                 gp = gpar(col = "#9CA3AF", lwd = 0.8))
    }
  }

  pushViewport(viewport(layout.pos.row = 7, x = 0.5, width = 0.96))
  grid.roundrect(gp = gpar(fill = "#E8F5EF", col = cols["Retained claim"], lwd = 1.2), r = unit(2.5, "mm"))
  grid.text("RETAINED CLAIM", x = 0.035, y = 0.76, just = "left",
            gp = gpar(fontfamily = "sans", fontsize = 7.2, fontface = "bold", col = cols["Retained claim"]))
  grid.text(wrap(stages$evidence[5], 59), x = 0.035, y = 0.43, just = "left",
            gp = gpar(fontfamily = "sans", fontsize = 8.2, fontface = "bold", col = "#111827", lineheight = 1.05))
  grid.text(wrap(stages$interpretation[5], 47), x = 0.965, y = 0.42, just = "right",
            gp = gpar(fontfamily = "sans", fontsize = 6.9, col = "#374151", lineheight = 1.08))
  popViewport()

  grid.text(
    "Conceptual synthesis; equal box dimensions do not encode effect size or evidential weight.",
    vp = viewport(layout.pos.row = 8), x = unit(0.02, "npc"), just = "left",
    gp = gpar(fontfamily = "sans", fontsize = 6.4, col = "#6B7280")
  )
  popViewport()
}

pdf(pdf_file, width = 183 / 25.4, height = 143 / 25.4, family = "sans", useDingbats = FALSE)
draw_figure()
dev.off()

png_type <- if (capabilities("aqua")) "quartz" else "cairo"
png(png_file, width = 183, height = 143, units = "mm", res = 600, type = png_type, bg = "white")
draw_figure()
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "Figure5_evidence_contraction_sessionInfo.txt"))
message("Wrote: ", png_file, ", ", pdf_file, ", and ", tsv_file)
