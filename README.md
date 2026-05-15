# CTEnrich

Cell type enrichment scoring and visualization for bulk RNA-seq.

Given a set of DESeq2 results (one per condition vs. a baseline) and a table of cell-type marker gene sets, CTEnrich scores each condition by averaging the log2 fold-change values of each cell type's markers. High mean LFC indicates the condition has upregulated that cell type's signature relative to the baseline.

> **Recommended baseline:** pluripotent cells (e.g. iPSCs or ESCs). Because pluripotent cells express few lineage-specific markers, virtually any differentiated signature will appear as a positive LFC, maximising the diversity of detectable enriched cell types. Publicly available bulk pluripotent RNA-seq data can be used if an in-house baseline is unavailable, but samples should be processed together in the same experiment to minimise batch effects.

## Installation

```r
# Development version
devtools::install_github("masayatakasaki/CTEnrich")

# Or load without installing (e.g. on HPC)
devtools::load_all("~/R/CTEnrich")
```

## Input format

### `res`
A DESeq2 `DESeqResults` object or any data frame with:
- `log2FoldChange` column
- Gene symbols as rownames

### `markers`
A data frame with:

| Column | Required | Description |
|---|---|---|
| `cell_type` | yes | Cell type name (must be unique) |
| `marker_genes` | yes | List-column of character vectors (gene symbols) |
| `germ_layer` | no | Used for colouring in plots |

Use `split_markers()` to convert a comma-separated string column to the required list-column format, and `validate_markers()` to check the input before scoring.

Gene symbols must be in the same namespace as `rownames(res)` — e.g. both HGNC human symbols, or both MGI mouse symbols.

## Built-in marker set

CTEnrich ships with `markers_human`, a set of 252 cell type marker gene lists
derived from the mouse developmental atlas of Qiu et al. (2024), converted to
human orthologs via Ensembl BioMart. It covers six germ layers (Ectoderm,
Endoderm, Extraembryonic, Mesoderm, Other, PGC) and is suitable for scoring
human bulk RNA-seq data against developmental cell type signatures.

```r
data(markers_human)
# 252 cell types, ~19 markers each
```

## Quick start

```r
library(CTEnrich)

# Use the built-in marker set (Qiu et al. 2024, human orthologs)
data(markers_human)

# Or supply your own — split_markers() converts comma-separated strings
markers <- data.frame(
  cell_type    = c("Cardiomyocyte", "Hepatocyte", "Neuron"),
  marker_genes = c("TNNT2, MYH7, ACTC1", "ALB, AFP, HNF4A", "MAP2, SYP, NCAM1"),
  germ_layer   = c("Mesoderm", "Endoderm", "Ectoderm")
)
markers <- split_markers(markers, "marker_genes", sep = ", ")
validate_markers(markers)

# Score all conditions (resList = named list of DESeq2 results)
markerResList <- lapply(resList, marker_enrich_table, markers = markers_human)

# Per-condition plots
gen_layer(markerResList[["condition_a"]], name = "condition_a", filter = TRUE)

# Cross-condition comparison
gen_comp_multi(markerResList, thresh = 2, figs_dir = "figs/")

# AUC progression across conditions
gen_auc_thresh(markerResList, thresh = 2, label = "AUC_comparison", figs_dir = "figs/")
gen_auc_bar(   markerResList, thresh = 2, label = "AUC_bar",        figs_dir = "figs/")
gen_n_enriched(markerResList, thresh = 2, label = "N_enriched",     figs_dir = "figs/")
```

## Function reference

| Function | Description |
|---|---|
| `marker_enrich_table(res, markers)` | Score all cell types for one condition |
| `validate_markers(markers)` | Check markers input format |
| `split_markers(df, col, sep)` | Convert comma-separated string column to list-column |
| `gen_layer(res, name, thresh, filter)` | Dot plot of mean LFC per cell type |
| `gen_comp_thresh(ordered_res_list, thresh)` | Compare conditions against an anchor |
| `gen_comp_multi(res_list, thresh)` | Run `gen_comp_thresh` for every condition as anchor |
| `compute_auc(res_list, thresh)` | Return named vector of AUC scores |
| `gen_auc_thresh(res_list, thresh, label)` | Smooth rank curves with AUC in caption |
| `gen_auc_bar(res_list, thresh, label)` | Bar plot of AUC per condition |
| `gen_n_enriched(res_list, thresh, label)` | Bar plot of enriched cell type count per condition |

All plot functions accept `figs_dir = NULL` to return the plot object without saving.

## References

Qiu C, et al. (2024) A single-cell time-lapse of mouse development from gastrulation to birth. *Nature*.
