# Signalling heatmaps across the five conditions, one page per condition.
#
# Each entry in SETS is either a single pathway or several summed together.
# netVisual_heatmap only accepts one pathway at a time, so a multi-pathway set
# is built by hand from netP$prob and styled to match.
#
#   Rscript cellchat_heatmaps.R

library(CellChat)
library(ComplexHeatmap)
library(circlize)

# EDIT THIS. R cannot read config.py, so the path is set here.
CELLCHAT_DIR <- "/path/to/CellChat/Nerandomilast_10Celltypes"

CONDS <- c("WT", "I73T_Week_4", "Nerandomilast",
           "I73T_Week_4_18M", "Nerandomilast_18M")

SETS <- list(
  BAFF           = "BAFF",
  TNFsuperfamily = c("TNF", "RANKL", "TWEAK", "BAFF", "LIGHT")
)

ccs <- lapply(CONDS, function(cn)
  readRDS(file.path(CELLCHAT_DIR, sprintf("cellchat_%s.rds", cn))))
names(ccs) <- CONDS


# Every figure is written twice, PDF for print and SVG for Illustrator.
to_pdf_and_svg <- function(base, width, height, draw_fn) {
  cairo_pdf(paste0(base, ".pdf"), width = width, height = height)
  draw_fn()
  dev.off()

  svg(paste0(base, ".svg"), width = width, height = height)
  draw_fn()
  dev.off()

  cat("saved:", paste0(base, ".pdf"), "and", paste0(base, ".svg"), "\n")
}


# --- one pathway: CellChat draws it -------------------------------------------
draw_single <- function(set_name, pathway) {
  function() {
    for (cn in CONDS) {
      cc <- ccs[[cn]]
      if (!(pathway %in% cc@netP$pathways)) {
        cat(sprintf("  %s: not detected\n", cn))
        next
      }
      cat(sprintf("  %s: %s\n", cn, pathway))
      draw(netVisual_heatmap(cc,
                             signaling     = pathway,
                             color.heatmap = "Reds",
                             title.name    = sprintf("%s - %s", set_name, cn)))
    }
  }
}


# --- several pathways: sum them, then draw ------------------------------------
# Each condition sums only the pathways it actually has, and the panel title
# records which those were. Colours match netVisual_heatmap - the same Reds
# ramp and the same per-cell-type bars down the side and across the top.
draw_aggregate <- function(set_name, pathways) {
  mats <- lapply(CONDS, function(cn) {
    cc <- ccs[[cn]]
    present <- intersect(pathways, cc@netP$pathways)
    cat(sprintf("  %s: %s\n", cn, paste(present, collapse = ", ")))
    if (length(present) == 0) return(NULL)
    list(mat = Reduce(`+`, lapply(present, function(p) cc@netP$prob[ , , p])),
         pathways = present)
  })
  names(mats) <- CONDS

  # one colour scale across every panel, so intensity means the same thing in
  # all of them - netVisual_heatmap scales each panel to itself, which is fine
  # for one condition but misleading side by side
  global_max <- max(sapply(mats, function(x)
    if (is.null(x)) 0 else max(x$mat, na.rm = TRUE)))
  col_fun <- colorRamp2(seq(0, global_max, length.out = 9),
                        RColorBrewer::brewer.pal(9, "Reds"))

  function() {
    for (cn in CONDS) {
      x <- mats[[cn]]
      if (is.null(x)) next
      m <- x$mat

      groups <- rownames(m)
      group_cols <- scPalette(length(groups))
      names(group_cols) <- groups

      hm <- Heatmap(
        m,
        name = "Comm. Prob.",
        col  = col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_side    = "left",
        column_names_side = "bottom",
        column_names_rot  = 45,
        row_names_gp      = gpar(fontsize = 8),
        column_names_gp   = gpar(fontsize = 8),
        column_title      = sprintf("%s - %s (%s)", set_name, cn,
                                    paste(x$pathways, collapse = "+")),
        column_title_gp   = gpar(fontsize = 11, fontface = "bold"),

        # rows send, columns receive
        top_annotation = HeatmapAnnotation(
          target = groups,
          incoming = anno_barplot(colSums(m), border = FALSE,
                                  gp = gpar(fill = group_cols, col = NA)),
          col = list(target = group_cols),
          show_legend = FALSE, show_annotation_name = FALSE,
          simple_anno_size = grid::unit(2, "mm")),

        left_annotation = rowAnnotation(
          source = groups,
          col = list(source = group_cols),
          show_legend = FALSE, show_annotation_name = FALSE,
          simple_anno_size = grid::unit(2, "mm")),

        right_annotation = rowAnnotation(
          outgoing = anno_barplot(rowSums(m), border = FALSE,
                                  gp = gpar(fill = group_cols, col = NA)),
          show_annotation_name = FALSE))

      draw(hm)
    }
  }
}


for (set_name in names(SETS)) {
  pathways <- SETS[[set_name]]
  cat(set_name, ":\n", sep = "")

  base <- file.path(CELLCHAT_DIR, sprintf("heatmap_%s_5conditions", set_name))

  draw_fn <- if (length(pathways) == 1) {
    draw_single(set_name, pathways)
  } else {
    draw_aggregate(set_name, pathways)
  }

  to_pdf_and_svg(base, width = 9, height = 8, draw_fn = draw_fn)
}
