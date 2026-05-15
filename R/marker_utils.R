# validate_markers: check that a markers data frame meets the input contract
# for marker_enrich_table.  Stops with a clear message on the first failure.
validate_markers <- function(markers) {
  if (!is.data.frame(markers))
    stop("`markers` must be a data frame", call. = FALSE)

  required <- c("cell_type", "marker_genes")
  missing  <- setdiff(required, colnames(markers))
  if (length(missing) > 0)
    stop("required columns missing from `markers`: ",
         paste(missing, collapse = ", "), call. = FALSE)

  if (!is.list(markers$marker_genes))
    stop("`markers$marker_genes` must be a list-column; ",
         "use split_markers() to convert a comma-separated string column",
         call. = FALSE)

  empty <- vapply(markers$marker_genes, function(g) length(g) == 0, logical(1))
  if (any(empty))
    warning(sum(empty), " cell type(s) have empty marker gene lists",
            call. = FALSE)

  dups <- markers$cell_type[duplicated(markers$cell_type)]
  if (length(dups) > 0)
    warning("duplicate cell types in `markers` — marker_enrich_table will keep ",
            "only the highest-scoring entry: ",
            paste(unique(dups), collapse = ", "), call. = FALSE)

  invisible(markers)
}

# split_markers: convert a comma-separated string column into the list-column
# format required by marker_enrich_table.
split_markers <- function(df, col, sep = ", ") {
  if (!col %in% colnames(df))
    stop("column '", col, "' not found in data frame", call. = FALSE)
  df[[col]] <- strsplit(as.character(df[[col]]), sep, fixed = TRUE)
  df
}
