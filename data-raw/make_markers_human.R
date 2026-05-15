# Generates markers_human.rda — the built-in CTEnrich marker dataset.
# Source: Qiu et al. 2024 mouse developmental atlas, markers converted to
# human orthologs via Ensembl BioMart (mart_export_human_mouse.txt).
#
# Run once to regenerate the package dataset:
#   Rscript data-raw/make_markers_human.R

EXCEL_FILE   <- "/net/hamazaki/vol1/project/masaya/SP8x8/data/markers/CellTypeEnrich.xlsx"
ORTHO_FILE   <- "/net/hamazaki/vol1/project/masaya/SP8x8/data/markers/mart_export_human_mouse.txt"
OUT_FILE     <- "data/markers_human.rda"

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# Load raw marker table (cell type, comma-separated mouse genes, germ layer)
marker_raw            <- as.data.frame(read_excel(EXCEL_FILE, sheet = "CXTable"))[, 2:4]
colnames(marker_raw)  <- c("cell_type", "marker_genes", "germ_layer")
marker_raw            <- distinct(marker_raw, cell_type, .keep_all = TRUE)

# Split comma-separated gene strings into character vectors
marker_raw$marker_genes <- strsplit(as.character(marker_raw$marker_genes), ", ", fixed = TRUE)

# Convert mouse gene symbols to human orthologs
link <- read.table(ORTHO_FILE, header = TRUE)
marker_raw$marker_genes <- lapply(marker_raw$marker_genes,
  function(gs) link$Human_name[link$Mouse_name %in% gs])

markers_human <- marker_raw

cat("Cell types:          ", nrow(markers_human), "\n")
cat("Avg markers per type:", round(mean(lengths(markers_human$marker_genes)), 1), "\n")
cat("Germ layers:         ", paste(sort(unique(markers_human$germ_layer)), collapse = ", "), "\n")

save(markers_human, file = OUT_FILE, compress = "xz")
cat("Saved to", OUT_FILE, "\n")
