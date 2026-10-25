library(grid)

out_dir <- file.path(getwd(), "CLAIM_CITATION_AUDIT", "aad_revision_candidate", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

draw_ga <- function() {
  grid.newpage()
  navy <- "#17365D"
  blue <- "#DCE8F5"
  green <- "#E3F0E8"
  sand <- "#F3EEE5"
  border <- "#64748B"
  accent <- "#2F6DA5"

  grid.text(
    "Cross-species evidence suggests a bounded arterial-remodelling overlap",
    x = 0.5, y = 0.93,
    gp = gpar(fontfamily = "Helvetica", fontsize = 21, fontface = "bold", col = navy)
  )

  xs <- c(0.18, 0.50, 0.82)
  fills <- c(blue, green, sand)
  heads <- c("Animal discovery", "Preliminary human support", "Boundary tests")
  body <- list(
    c("549 shared genes", "478 concordant; 71 discordant", "155 shared pathways"),
    c("Ossification: 3 of 3 datasets", "Matrix and smooth-muscle signals varied", "5, 10 and 15 cases across datasets"),
    c("No renal trait passed BH correction", "Genetic correlation was non-estimable", "Classifier signal lacked specificity")
  )

  for (i in seq_along(xs)) {
    grid.roundrect(
      x = xs[i], y = 0.58, width = 0.27, height = 0.43, r = unit(4, "mm"),
      gp = gpar(fill = fills[i], col = border, lwd = 1.2)
    )
    grid.text(heads[i], x = xs[i], y = 0.72,
              gp = gpar(fontfamily = "Helvetica", fontsize = 15, fontface = "bold", col = navy))
    yy <- c(0.61, 0.555, 0.50)
    for (j in 1:3) {
      grid.text(body[[i]][j], x = xs[i], y = yy[j],
                gp = gpar(fontfamily = "Helvetica", fontsize = 10.8, col = "#222222"))
    }
  }

  grid.text("Uraemic calcification and\ndissection-or-aneurysm models",
            x = xs[1], y = 0.415,
            gp = gpar(fontfamily = "Helvetica", fontsize = 9, col = "#222222"))

  arrow_gp <- gpar(col = accent, fill = accent, lwd = 2)
  grid.lines(x = unit(c(0.325, 0.365), "npc"), y = unit(c(0.58, 0.58), "npc"),
             arrow = arrow(type = "closed", length = unit(2.6, "mm")), gp = arrow_gp)
  grid.lines(x = unit(c(0.645, 0.685), "npc"), y = unit(c(0.58, 0.58), "npc"),
             arrow = arrow(type = "closed", length = unit(2.6, "mm")), gp = arrow_gp)

  grid.roundrect(x = 0.5, y = 0.20, width = 0.73, height = 0.17,
                 r = unit(4, "mm"), gp = gpar(fill = navy, col = navy))
  grid.text("Most reproducible human signal: ossification-related expression",
            x = 0.5, y = 0.235,
            gp = gpar(fontfamily = "Helvetica", fontsize = 14.5, fontface = "bold", col = "white"))
  grid.text("The overlap guides focused studies but does not establish causality or nominate treatment targets",
            x = 0.5, y = 0.165,
            gp = gpar(fontfamily = "Helvetica", fontsize = 10.8, col = "white"))
}

stopifnot(capabilities("aqua"))
grDevices::quartz(type="pdf",file=file.path(out_dir, "GraphicalAbstract_R1.pdf"), width = 13, height = 5,
    family = "Helvetica")
draw_ga(); dev.off()

# Render both derivatives from the new PDF; fail rather than retain stale images.
pdf_path <- file.path(out_dir, "GraphicalAbstract_R1.pdf")
stopifnot(system2("/opt/homebrew/bin/pdftocairo", c("-svg", shQuote(pdf_path),
  shQuote(file.path(out_dir,"GraphicalAbstract_R1.svg")))) == 0)
stopifnot(system2("/opt/homebrew/bin/pdftoppm", c("-png", "-singlefile", "-r", "600",
  shQuote(pdf_path), shQuote(file.path(out_dir,"GraphicalAbstract_R1")))) == 0)
