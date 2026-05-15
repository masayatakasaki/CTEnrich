# Score cell type enrichment and visualize results across conditions.
#
# Input contract for `markers`:
#   cell_type    – character, must be unique (duplicates resolved by keeping
#                  the highest-scoring row after scoring)
#   marker_genes – list-column of character vectors; gene symbols must be in
#                  the same namespace as rownames(res) (e.g. both human)
#   germ_layer   – character, optional; used for colour in layer plots
#   any other columns pass through to the output table unchanged
#
# Use validate_markers() and split_markers() in marker_utils.R to prepare input.

# For each cell type, average the LFC values of its marker genes found in res.
# n = genes matched; n_total = genes in the marker set (coverage = n/n_total).
# Rows are ordered by mean descending; duplicated cell types are collapsed to
# the highest-scoring entry (a side effect of input duplicates, not intentional
# deduplication — validate_markers() warns if this will happen).
marker_enrich_table <- function(res, markers) {
  res_df <- as.data.frame(res)
  if (!"log2FoldChange" %in% colnames(res_df))
    stop("`res` must have a 'log2FoldChange' column", call. = FALSE)
  if (is.null(rownames(res_df)) || all(rownames(res_df) == as.character(seq_len(nrow(res_df)))))
    stop("`res` must have gene symbols as rownames", call. = FALSE)
  if (!is.list(markers$marker_genes))
    stop("`markers$marker_genes` must be a list-column of character vectors", call. = FALSE)

  lfc <- setNames(res_df$log2FoldChange, rownames(res_df))

  stats <- lapply(markers$marker_genes, function(gs) {
    vals <- lfc[intersect(gs, names(lfc))]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0)
      warning("a cell type has 0 marker genes in the result — mean will be NaN",
              call. = FALSE)
    list(mean = mean(vals), median = median(vals), n = length(vals),
         n_total = length(gs), stdev = sd(vals), enrich_vals = unname(vals))
  })

  out <- dplyr::bind_cols(markers, dplyr::tibble(
    mean        = vapply(stats, `[[`, numeric(1),   "mean"),
    median      = vapply(stats, `[[`, numeric(1),   "median"),
    n           = vapply(stats, `[[`, integer(1),   "n"),
    n_total     = vapply(stats, `[[`, integer(1),   "n_total"),
    stdev       = vapply(stats, `[[`, numeric(1),   "stdev"),
    enrich_vals = lapply(stats, `[[`, "enrich_vals")
  ))
  out <- out[order(out$mean, decreasing = TRUE), ]
  dplyr::distinct(out, cell_type, .keep_all = TRUE)
}

# Reorder rows of df to match ref_levels order.
# Unmatched ref_levels produce NA rows — intentional so plot_comp_thresh can
# plot all base cell types even when a comparison condition is missing some.
reorder_like <- function(df, ref_levels) {
  idx     <- match(ref_levels, as.character(df$cell_type))
  missing <- ref_levels[is.na(idx)]
  if (length(missing) > 0)
    warning("cell types in ref_levels not found in df: ",
            paste(missing, collapse = ", "), call. = FALSE)
  df[idx, , drop = FALSE]
}

# Produce all cyclic rotations of a named list so that each element appears
# first exactly once — used by plot_comp_multi to cycle every condition as base.
rotate_list <- function(x) {
  n <- length(x)
  lapply(seq_len(n), function(i) x[c(seq(i, n), if (i > 1) seq_len(i - 1))])
}

# Dot plot of mean LFC per cell type with individual marker gene LFCs jittered.
# filter = FALSE: all cell types (overview); filter = TRUE: above thresh only.
# thresh = NA with filter = TRUE shows the top 30 by mean instead.
# germ_layer column is optional; if absent all points use the default colour.
# coord_flip(ylim=...) is used instead of ylim() to avoid dropping data before
# geom_jitter runs (ylim() removes points, coord clipping does not).
plot_layer <- function(res, name, thresh = 2, filter = TRUE, figs_dir = ".") {
  res_f <- if (!filter) {
    res
  } else if (is.na(thresh)) {
    head(res, 30)
  } else {
    res[!is.na(res$mean) & res$mean > thresh, ]
  }

  if (nrow(res_f) == 0) {
    message("No cell types above threshold for: ", name, " — skipping")
    return(invisible(NULL))
  }

  germ_col <- if ("germ_layer" %in% colnames(res_f)) "germ_layer" else NULL

  unnest_df <- res_f |>
    dplyr::select(dplyr::any_of(c("cell_type", "germ_layer", "enrich_vals"))) |>
    tidyr::unnest_longer(enrich_vals)

  if (!is.null(germ_col))
    unnest_df <- dplyr::mutate(unnest_df,
                               germ_layer = dplyr::coalesce(germ_layer, "Other"))

  p <- res_f |>
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type))))

  if (!is.null(germ_col))
    p <- dplyr::mutate(p, germ_layer = dplyr::coalesce(germ_layer, "Other"))

  label_size <- if (filter) max(4, min(8, 120 / nrow(res_f))) else NULL

  p <- ggplot2::ggplot(p, ggplot2::aes(cell_type, mean)) +
    ggplot2::geom_jitter(
      data = unnest_df,
      ggplot2::aes(y = enrich_vals,
                   color = if (!is.null(germ_col)) germ_layer else NULL),
      width = 0.15, alpha = 0.5, size = 1, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = if (!is.null(germ_col)) germ_layer else NULL), size = 3
    ) +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed") +
    ggplot2::coord_flip(ylim = c(-10, 10)) +
    ggplot2::theme_bw() +
    if (!is.null(label_size))
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = label_size))

  caption <- if (filter && !is.na(thresh))
    paste0(nrow(res_f), " cell types, threshold = ", thresh) else NULL

  p <- p + ggplot2::labs(
    title   = paste0(name, if (filter) " enriched cell types" else " all cell types"),
    caption = caption,
    y = "mean log2FC", color = NULL
  )

  if (!is.null(figs_dir)) {
    suffix <- if (filter) "LayerThresh" else "LayerFull"
    h      <- if (filter) 6 else 25
    ggplot2::ggsave(file.path(figs_dir, paste0(name, suffix, ".tiff")), p, height = h)
  }
  invisible(p)
}

# Cross-condition comparison dot plot. The first element of ordered_res_list is
# the anchor — its cell types define the row set and order. All other conditions
# are projected onto those rows via reorder_like(). The anchor is placed last so
# its points render on top.
# filter = TRUE: anchor filtered to above-thresh cell types; anchor highlighted
#   with a binary colour; height = 6; saves as {name}CompThresh.tiff.
# filter = FALSE: all anchor cell types shown; each condition a distinct colour;
#   height = 25; saves as {name}CompFull.tiff.
plot_comp <- function(ordered_res_list, thresh = 2, filter = TRUE, figs_dir = ".") {
  base_name <- names(ordered_res_list)[1]
  base_df   <- ordered_res_list[[1]]

  if (filter) {
    base_df <- base_df[!is.na(base_df$mean) & base_df$mean > thresh, ]
    if (nrow(base_df) == 0) {
      message("'", base_name, "' has 0 cell types above threshold — skipping")
      return(invisible(NULL))
    }
  }

  clusters  <- as.character(base_df$cell_type)
  all_list  <- c(lapply(ordered_res_list[-1], reorder_like, clusters),
                 ordered_res_list[1])

  plot_data <- dplyr::bind_rows(lapply(names(all_list), function(nm) {
    data.frame(cell_type = as.character(all_list[[nm]]$cell_type),
               condition = nm,
               mean      = all_list[[nm]]$mean)
  }))

  label_size <- max(4, min(8, 120 / length(clusters)))

  p <- plot_data |>
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(clusters))) |>
    ggplot2::ggplot(ggplot2::aes(
      cell_type, mean,
      color = if (filter) condition == base_name else condition
    )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed") +
    ggplot2::coord_flip(ylim = if (filter) c(-3, 10) else c(-10, 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = label_size)) +
    ggplot2::labs(
      title   = paste0(base_name, if (filter) " comparison" else " all cell types comparison"),
      caption = paste0("threshold: ", thresh),
      y = "mean log2FC"
    )

  if (!is.null(figs_dir)) {
    suffix <- if (filter) "CompThresh" else "CompFull"
    h      <- if (filter) 6 else 25
    ggplot2::ggsave(file.path(figs_dir, paste0(base_name, suffix, ".tiff")),
                    p, height = h, width = 12)
  }
  invisible(p)
}

# Run plot_comp with each condition as the anchor in turn.
plot_comp_multi <- function(res_list, thresh = 2, filter = TRUE, figs_dir = ".") {
  invisible(lapply(rotate_list(res_list), plot_comp,
                   thresh = thresh, filter = filter, figs_dir = figs_dir))
}

# Uniqueness score per cell type = anchor_mean - max(all_other_means).
# Positive values indicate cell types whose signature is stronger in the anchor
# than in any comparison condition. Top `top` cell types shown.
# Requires at least 2 conditions; cell types with no data in any other condition
# are excluded (NA uniqueness).
plot_uniqueness <- function(ordered_res_list, top = 30, figs_dir = ".") {
  if (length(ordered_res_list) < 2)
    stop("`ordered_res_list` must have at least 2 conditions", call. = FALSE)

  base_name <- names(ordered_res_list)[1]
  base_df   <- ordered_res_list[[1]]
  clusters  <- as.character(base_df$cell_type)

  other_means <- do.call(cbind, lapply(
    lapply(ordered_res_list[-1], reorder_like, clusters), `[[`, "mean"
  ))
  max_other <- apply(other_means, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
  })

  germ_col  <- if ("germ_layer" %in% colnames(base_df)) "germ_layer" else NULL
  plot_data <- data.frame(cell_type = clusters, us = base_df$mean - max_other,
                          stringsAsFactors = FALSE)
  if (!is.null(germ_col))
    plot_data$germ_layer <- dplyr::coalesce(base_df$germ_layer, "Other")

  plot_data <- plot_data[!is.na(plot_data$us), ]
  plot_data <- head(plot_data[order(plot_data$us, decreasing = TRUE), ], top)

  if (nrow(plot_data) == 0) {
    message("No uniqueness scores computable for: ", base_name, " — skipping")
    return(invisible(NULL))
  }

  p <- plot_data |>
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) |>
    ggplot2::ggplot(ggplot2::aes(cell_type, us,
                                  color = if (!is.null(germ_col)) germ_layer else NULL)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0(base_name, " uniqueness scores"),
                  y = "uniqueness (log2FC)", color = NULL)

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, paste0(base_name, "Uniqueness.tiff")),
                    p, height = 6)
  invisible(p)
}

# Run plot_uniqueness with each condition as the anchor in turn.
plot_uniqueness_multi <- function(res_list, top = 30, figs_dir = ".") {
  invisible(lapply(rotate_list(res_list), plot_uniqueness,
                   top = top, figs_dir = figs_dir))
}

# Bar chart of mean developmental time of enriched cell types per condition.
# Requires a 'time' column in the markers data frame (passed through by
# marker_enrich_table). Conditions with no enriched cell types that have time
# annotations are silently dropped.
plot_mean_time <- function(res_list, thresh = 2, label, figs_dir = ".") {
  has_time <- vapply(res_list, function(r) "time" %in% colnames(r), logical(1))
  if (!all(has_time))
    stop("All elements of `res_list` must have a 'time' column. ",
         "Add a 'time' column to your markers data frame before calling ",
         "marker_enrich_table().", call. = FALSE)

  mean_times <- vapply(names(res_list), function(nm) {
    r     <- res_list[[nm]]
    above <- r[!is.na(r$mean) & r$mean >= thresh & !is.na(r$time), ]
    if (nrow(above) == 0) NA_real_ else mean(above$time)
  }, numeric(1))

  plot_data <- data.frame(condition = names(mean_times), mean_time = mean_times,
                          row.names = NULL)
  plot_data <- plot_data[!is.na(plot_data$mean_time), ]

  if (nrow(plot_data) == 0) {
    message("No conditions have enriched cell types with time annotations — skipping")
    return(invisible(NULL))
  }

  plot_data <- plot_data[order(plot_data$mean_time), ]
  plot_data$condition <- factor(plot_data$condition, levels = plot_data$condition)

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(condition, mean_time, fill = condition)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::labs(y       = "mean embryonic day of enriched cell types",
                  caption = paste0("threshold: ", thresh))

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, paste0(label, "MeanTime.tiff")), p, width = 8)
  invisible(p)
}

# For consecutively-ordered conditions, plots cell types newly enriched at each
# step (above thresh in condition i but not i-1) by developmental time vs
# log2(mean LFC). Requires a 'time' column in the markers data frame and the
# ggrepel package for labels.
plot_order_comp <- function(res_list, thresh = 2, figs_dir = ".") {
  has_time <- vapply(res_list, function(r) "time" %in% colnames(r), logical(1))
  if (!all(has_time))
    stop("All elements of `res_list` must have a 'time' column. ",
         "Add a 'time' column to your markers data frame before calling ",
         "marker_enrich_table().", call. = FALSE)
  if (length(res_list) < 2)
    stop("`res_list` must have at least 2 conditions", call. = FALSE)

  nms       <- names(res_list)
  plot_data <- dplyr::bind_rows(lapply(seq(2, length(nms)), function(i) {
    before_types <- as.character(
      res_list[[nms[i - 1]]]$cell_type[
        !is.na(res_list[[nms[i - 1]]]$mean) &
          res_list[[nms[i - 1]]]$mean >= thresh
      ]
    )
    after <- res_list[[nms[i]]]
    after <- after[!is.na(after$mean) & after$mean >= thresh, ]
    newly <- after[!as.character(after$cell_type) %in% before_types, ]
    if (nrow(newly) == 0) return(NULL)

    out <- data.frame(cell_type = as.character(newly$cell_type),
                      condition = nms[i],
                      mean      = newly$mean,
                      time      = newly$time,
                      stringsAsFactors = FALSE)
    if ("germ_layer" %in% colnames(newly)) out$germ_layer <- newly$germ_layer
    out
  }))

  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("No newly enriched cell types found across conditions — skipping")
    return(invisible(NULL))
  }

  germ_col <- if ("germ_layer" %in% colnames(plot_data)) "germ_layer" else NULL

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(time, log2(mean), color = condition,
                 shape = if (!is.null(germ_col)) germ_layer else NULL)
  ) +
    ggplot2::geom_point(size = 4) +
    ggrepel::geom_text_repel(ggplot2::aes(label = cell_type),
                              max.overlaps = Inf, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "mouse embryonic day", y = "log2 mean LFC",
                  shape = NULL, color = NULL)

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, "OrderComp.tiff"), p, width = 10)
  invisible(p)
}

# Raw sum of mean LFC values >= thresh per condition.
# NAs (cell types absent from res) are excluded, not treated as 0.
compute_auc <- function(res_list, thresh = 2) {
  vapply(names(res_list), function(nm) {
    m <- res_list[[nm]]$mean
    sum(m[!is.na(m) & m >= thresh])
  }, numeric(1))
}

# Smooth curves of mean LFC vs cell-type rank, one per condition.
# The LOESS curve and the AUC number are intentionally separate: the curve
# shows the trend shape; the AUC (raw sum, in caption) measures total area.
plot_auc_curve <- function(res_list, thresh = 2, label, figs_dir = ".") {
  auc_vals <- compute_auc(res_list, thresh)

  auc_data <- dplyr::bind_rows(lapply(names(res_list), function(nm) {
    m <- res_list[[nm]]$mean
    m[is.na(m) | m < thresh] <- 0
    data.frame(cell_type = as.character(res_list[[nm]]$cell_type),
               condition = nm,
               mean      = m,
               rank      = seq_len(nrow(res_list[[nm]])))
  }))

  auc_caption <- paste(paste0(names(auc_vals), ": ", round(auc_vals, 1)),
                       collapse = "\n")

  p <- ggplot2::ggplot(auc_data, ggplot2::aes(rank, mean, color = condition)) +
    ggplot2::geom_smooth(se = FALSE, ggplot2::aes(linetype = condition)) +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position        = "inside",
                   legend.position.inside = c(.75, .85),
                   legend.text            = ggplot2::element_text(size = 8),
                   legend.title           = ggplot2::element_blank()) +
    ggplot2::labs(x = "cell type rank", y = "mean log2FC",
                  caption = paste0("AUC (raw sum above thresh):\n", auc_caption)) +
    ggplot2::coord_cartesian(ylim = c(thresh, 7))

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, paste0(label, ".tiff")), p, width = 8)
  invisible(p)
}

# Bar plot of number of cell types with mean >= thresh per condition.
# Complements plot_auc_bar: breadth (how many) vs. magnitude (how strongly).
plot_n_enriched <- function(res_list, thresh = 2, label, figs_dir = ".") {
  counts <- vapply(names(res_list), function(nm) {
    m <- res_list[[nm]]$mean
    sum(!is.na(m) & m >= thresh)
  }, integer(1))

  plot_data <- data.frame(
    condition = factor(names(counts), levels = names(counts)),
    n         = counts
  )

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(condition, n, fill = condition)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::labs(y = "number of enriched cell types",
                  caption = paste0("threshold: ", thresh))

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, paste0(label, ".tiff")), p, width = 8)
  invisible(p)
}

# Bar plot of AUC score per condition (total enrichment magnitude).
# Complements plot_n_enriched: magnitude (how strongly) vs. breadth (how many).
plot_auc_bar <- function(res_list, thresh = 2, label, figs_dir = ".") {
  auc_vals <- compute_auc(res_list, thresh)
  plot_data <- data.frame(
    condition = factor(names(auc_vals), levels = names(auc_vals)),
    auc       = auc_vals
  )

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(condition, auc, fill = condition)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::labs(y = "cell type diversity (AUC)",
                  caption = paste0("threshold: ", thresh))

  if (!is.null(figs_dir))
    ggplot2::ggsave(file.path(figs_dir, paste0(label, ".tiff")), p, width = 8)
  invisible(p)
}
