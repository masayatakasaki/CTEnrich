# Generates markers_mouse.rda and markers_human.rda — the built-in CTEnrich
# marker datasets. Source: Qiu et al. 2024 mouse developmental atlas.
# Human orthologs mapped via Ensembl BioMart (mart_export_human_mouse.txt).
# Developmental time annotations from typetime_df.rds.
#
# Run once to regenerate (from package root):
#   Rscript data-raw/make_markers.R

EXCEL_FILE    <- "/net/hamazaki/vol1/project/masaya/SP8x8/data/markers/CellTypeEnrich.xlsx"
ORTHO_FILE    <- "~/R/bulk_refs/mart_export_human_mouse.txt"
TYPETIME_FILE <- "~/R/bulk_refs/typetime_df.rds"
OUT_MOUSE     <- "data/markers_mouse.rda"
OUT_HUMAN     <- "data/markers_human.rda"

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# Load raw marker table (cell type, comma-separated mouse genes, germ layer)
marker_raw           <- as.data.frame(read_excel(EXCEL_FILE, sheet = "CXTable"))[, 2:4]
colnames(marker_raw) <- c("cell_type", "marker_genes", "germ_layer")
marker_raw           <- distinct(marker_raw, cell_type, .keep_all = TRUE)
marker_raw$marker_genes <- strsplit(as.character(marker_raw$marker_genes), ", ", fixed = TRUE)

# Add developmental time annotations
typetime        <- readRDS(path.expand(TYPETIME_FILE))
marker_raw$time <- typetime$time[match(marker_raw$cell_type, typetime$celltype)]

# Mouse version — original MGI symbols
markers_mouse <- marker_raw

cat("--- markers_mouse ---\n")
cat("Cell types:          ", nrow(markers_mouse), "\n")
cat("Avg markers per type:", round(mean(lengths(markers_mouse$marker_genes)), 1), "\n")
cat("Germ layers:         ", paste(sort(unique(markers_mouse$germ_layer)), collapse = ", "), "\n")
cat("Time range:          ", range(markers_mouse$time, na.rm = TRUE), "\n")
save(markers_mouse, file = OUT_MOUSE, compress = "xz")
cat("Saved to", OUT_MOUSE, "\n\n")

# Human version — mouse symbols converted to human orthologs via BioMart
link <- read.table(path.expand(ORTHO_FILE), header = TRUE)
marker_raw$marker_genes <- lapply(marker_raw$marker_genes,
  function(gs) link$Human_name[link$Mouse_name %in% gs])

markers_human <- marker_raw

cat("--- markers_human ---\n")
cat("Cell types:          ", nrow(markers_human), "\n")
cat("Avg markers per type:", round(mean(lengths(markers_human$marker_genes)), 1), "\n")
cat("Germ layers:         ", paste(sort(unique(markers_human$germ_layer)), collapse = ", "), "\n")
cat("Time range:          ", range(markers_human$time, na.rm = TRUE), "\n")
save(markers_human, file = OUT_HUMAN, compress = "xz")
cat("Saved to", OUT_HUMAN, "\n")
