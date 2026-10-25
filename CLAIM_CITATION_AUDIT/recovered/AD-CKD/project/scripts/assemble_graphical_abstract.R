#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: assemble_graphical_abstract.R <output-directory>")
out <- normalizePath(args[[1]], mustWork = TRUE)

draw_ga <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("Cross-species evidence narrows to a bounded arterial-remodelling state",
                  x = 0.5, y = 0.91,
                  gp = grid::gpar(fontfamily = "sans", fontsize = 16, fontface = "bold", col = "#17365D"))

  boxes <- data.frame(
    x = c(0.18, 0.50, 0.82),
    title = c("Animal discovery", "Human arterial support", "Boundary tests"),
    body = c(
      "549 shared genes\n155 shared pathways\nMatrix + smooth-muscle programmes",
      "Narrower, context-dependent signal\nMatrix organization + ossification\nSame-site dataset: n = 5 + 5",
      "No renal trait passed BH correction\nGenetic evidence remained limited\nImmune/classifier signals did not generalize"
    ),
    fill = c("#DCEAF7", "#E6F1E8", "#F3EEE4"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(boxes))) {
    grid::grid.roundrect(x = boxes$x[i], y = 0.57, width = 0.25, height = 0.43,
                         r = grid::unit(3, "mm"),
                         gp = grid::gpar(fill = boxes$fill[i], col = "#6A7785", lwd = 1.2))
    grid::grid.text(boxes$title[i], x = boxes$x[i], y = 0.71,
                    gp = grid::gpar(fontfamily = "sans", fontsize = 12, fontface = "bold", col = "#17365D"))
    grid::grid.text(boxes$body[i], x = boxes$x[i], y = 0.54,
                    gp = grid::gpar(fontfamily = "sans", fontsize = 9.5, col = "#202020", lineheight = 1.25))
  }
  for (x in c(0.325, 0.645)) {
    grid::grid.lines(x = c(x, x + 0.03), y = c(0.57, 0.57),
                     arrow = grid::arrow(type = "closed", length = grid::unit(3, "mm")),
                     gp = grid::gpar(col = "#4472A4", lwd = 2))
  }
  grid::grid.roundrect(x = 0.5, y = 0.18, width = 0.72, height = 0.17,
                       r = grid::unit(3, "mm"),
                       gp = grid::gpar(fill = "#17365D", col = NA))
  grid::grid.text("Shared arterial-wall molecular states guide focused mechanism studies",
                  x = 0.5, y = 0.205,
                  gp = grid::gpar(fontfamily = "sans", fontsize = 11.5, fontface = "bold", col = "white"))
  grid::grid.text("They do not establish a renal-to-dissection mechanism or nominate treatment targets",
                  x = 0.5, y = 0.145,
                  gp = grid::gpar(fontfamily = "sans", fontsize = 9.5, col = "white"))
}

grDevices::pdf(file.path(out, "AD_CKD_graphical_abstract.pdf"), width = 13 / 2.54, height = 5 / 2.54, family = "Helvetica", useDingbats = FALSE)
draw_ga(); grDevices::dev.off()

grDevices::png(file.path(out, "AD_CKD_graphical_abstract.png"), width = 2600, height = 1000, res = 200, type = "quartz")
draw_ga(); grDevices::dev.off()

grDevices::tiff(file.path(out, "AD_CKD_graphical_abstract.tiff"), width = 2600, height = 1000, res = 200, compression = "lzw", type = "quartz")
draw_ga(); grDevices::dev.off()
