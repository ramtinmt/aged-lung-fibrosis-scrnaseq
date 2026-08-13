# Circle and heatmap figures for signalling pathways across the five
# conditions, from the objects cellchat_run.R writes.
#
# Each entry in SETS is either a single pathway or several summed together.
# For every set this writes:
#
#   circle_<set>_5conditions.pdf      five panels side by side
#   heatmap_<set>_5conditions.pdf     one panel per condition
#
#   Rscript cellchat_figures.R

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(circlize)
})

# EDIT THIS. R cannot read config.py, so the path is set here.
CELLCHAT_DIR <- "/path/to/CellChat/Mes_Immune_FDC"

CONDS <- c("WT", "I73T_Week_4", "Nerandomilast",
           "I73T_Week_4_18M", "Nerandomilast_18M")

SETS <- list(
  TNF            = "TNF",
  RANKL          = "RANKL",
  TWEAK          = "TWEAK",
  BAFF           = "BAFF",
  LIGHT          = "LIGHT",
  TNFsuperfamily = c("TNF", "RANKL", "TWEAK", "BAFF", "LIGHT")
)

ccs <- lapply(CONDS, function(cn)
  readRDS(file.path(CELLCHAT_DIR, sprintf("cellchat_%s.rds", cn))))
names(ccs) <- CONDS


# Summed communication probability per condition. Not every member of a set is
# detected in every condition, so each sums only what it has - and the panel
# title records which, so a difference between panels can be checked against
# what went into them.
mats_for <- function(pathways) {
  out <- lapply(CONDS, function(cn) {
    cc <- ccs[[cn]]
    present <- intersect(pathways, cc@netP$pathways)
    if (length(present) == 0) return(NULL)
    list(mat = Reduce(`+`, lapply(present, function(p) cc@netP$prob[ , , p])),
         pathways = present)
  })
  names(out) <- CONDS
  out
}

to_pdf <- function(path, width, height, draw_fn) {
  cairo_pdf(path, width = width, height = height)
  draw_fn(); dev.off()
  cat("  ", basename(path), "\n", sep = "")
}


for (set_name in names(SETS)) {
  pathways <- SETS[[set_name]]
  single   <- length(pathways) == 1
  mats     <- mats_for(pathways)

  cat("\n", set_name, "\n", sep = "")
  for (cn in CONDS) {
    x <- mats[[cn]]
    cat(sprintf("  %-20s %s\n", cn,
                if (is.null(x)) "not detected" else paste(x$pathways, collapse = "+")))
  }

  # one scale per set, shared by its five panels: thickness and colour mean the
  # same thing in every condition. Sharing across sets instead would make the
  # weaker pathways invisible next to TNF.
  gmax <- max(sapply(mats, function(x)
    if (is.null(x)) 0 else max(x$mat, na.rm = TRUE)))

  # ---- circle ----
  to_pdf(file.path(CELLCHAT_DIR, sprintf("circle_%s_5conditions.pdf", set_name)),
         width = 30, height = 8, draw_fn = function() {
    par(mfrow = c(1, length(CONDS)), xpd = TRUE)
    for (cn in CONDS) {
      x <- mats[[cn]]
      if (is.null(x)) {
        plot.new(); title(sprintf("%s - %s\n(not detected)", set_name, cn))
        next
      }
      if (single) {
        netVisual_aggregate(ccs[[cn]],
                            signaling       = pathways,
                            layout          = "circle",
                            edge.weight.max = gmax,
                            signaling.name  = sprintf("%s - %s", set_name, cn))
      } else {
        netVisual_circle(x$mat,
                         vertex.weight   = as.numeric(table(ccs[[cn]]@idents)),
                         weight.scale    = TRUE,
                         edge.weight.max = gmax,
                         title.name      = sprintf("%s - %s", set_name, cn))
      }
    }
  })

  # ---- heatmap ----
  # netVisual_heatmap takes one pathway at a time, so a summed set has to be
  # drawn from the matrix directly
  col_fun <- colorRamp2(seq(0, gmax, length.out = 9),
                        RColorBrewer::brewer.pal(9, "Reds"))

  to_pdf(file.path(CELLCHAT_DIR, sprintf("heatmap_%s_5conditions.pdf", set_name)),
         width = 9, height = 8, draw_fn = function() {
    for (cn in CONDS) {
      x <- mats[[cn]]
      if (is.null(x)) next

      if (single) {
        draw(netVisual_heatmap(ccs[[cn]],
                               signaling     = pathways,
                               color.heatmap = "Reds",
                               title.name    = sprintf("%s - %s", set_name, cn)))
      } else {
        m <- x$mat
        groups <- rownames(m)
        group_cols <- scPalette(length(groups))
        names(group_cols) <- groups

        draw(Heatmap(
          m, col = col_fun, name = "Comm. Prob.",
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_names_side = "left", column_names_side = "bottom",
          column_names_rot = 45,
          row_names_gp    = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          column_title    = sprintf("%s - %s (%s)", set_name, cn,
                                    paste(x$pathways, collapse = "+")),
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),

          # rows send, columns receive
          top_annotation = HeatmapAnnotation(
            target = groups,
            incoming = anno_barplot(colSums(m), border = FALSE,
                                    gp = gpar(fill = group_cols, col = NA)),
            col = list(target = group_cols),
            show_legend = FALSE, show_annotation_name = FALSE,
            simple_anno_size = grid::unit(2, "mm")),
          left_annotation = rowAnnotation(
            source = groups, col = list(source = group_cols),
            show_legend = FALSE, show_annotation_name = FALSE,
            simple_anno_size = grid::unit(2, "mm")),
          right_annotation = rowAnnotation(
            outgoing = anno_barplot(rowSums(m), border = FALSE,
                                    gp = gpar(fill = group_cols, col = NA)),
            show_annotation_name = FALSE)))
      }
    }
  })
}

cat("\nfigures written to", CELLCHAT_DIR, "\n")
