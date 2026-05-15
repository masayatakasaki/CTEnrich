# CTEnrich

Cell type enrichment scoring and visualization for bulk RNA-seq.

Given a set of DESeq2 results (one per condition vs. a baseline) and a table of cell-type marker gene sets, CTEnrich scores each condition by averaging the log2 fold-change values of each cell type's markers. High mean LFC indicates the condition has upregulated that cell type's signature relative to the baseline.

> **Recommended baseline:** pluripotent cells (iPSCs or ESCs), which maximise detectable cell type diversity. Same-experiment, same-library-prep samples are ideal; public transcriptomes can be used in a pinch but may introduce batch effects.

## Installation

```r
devtools::install_github("masayatakasaki/CTEnrich")
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
| `time` | no | Developmental time; required for `plot_mean_time` and `plot_order_comp` |

Marker sets should ideally come from scRNA-seq data relevant to your biological context (e.g. top DE genes per cluster from an atlas). Context-matched markers will improve sensitivity. A built-in set derived from a mouse developmental atlas (converted to human orthologs) is included; see [Built-in marker set](#built-in-marker-set).

Use `split_markers()` to convert a comma-separated string column to the required list-column format, and `validate_markers()` to check the input before scoring. Any extra columns pass through `marker_enrich_table()` unchanged.

Gene symbols must be in the same namespace as `rownames(res)` — e.g. both HGNC human symbols, or both MGI mouse symbols.

## Built-in marker sets

CTEnrich ships with two marker sets derived from the mouse developmental atlas of Qiu et al. (2024), covering 252 cell types across six germ layers (Ectoderm, Endoderm, Extraembryonic, Mesoderm, Other, PGC) with ~19 markers each.

| Dataset | Symbols | Use with |
|---|---|---|
| `markers_human` | HGNC (human orthologs via Ensembl BioMart) | Human bulk RNA-seq |
| `markers_mouse` | MGI (original mouse symbols) | Mouse bulk RNA-seq |

```r
data(markers_human)  # for human data
data(markers_mouse)  # for mouse data
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
plot_layer(markerResList[["condition_a"]], name = "condition_a", filter = TRUE)

# Cross-condition comparison (filter = TRUE: above-threshold only; FALSE: all cell types)
plot_comp_multi(markerResList, thresh = 2, filter = TRUE,  figs_dir = "figs/")
plot_comp_multi(markerResList, thresh = 2, filter = FALSE, figs_dir = "figs/")

# Uniqueness scores (anchor mean minus max of all other conditions)
plot_uniqueness_multi(markerResList, top = 30, figs_dir = "figs/")

# AUC and summary plots
plot_auc_curve( markerResList, thresh = 2, label = "AUC_comparison", figs_dir = "figs/")
plot_auc_bar(   markerResList, thresh = 2, label = "AUC_bar",        figs_dir = "figs/")
plot_n_enriched(markerResList, thresh = 2, label = "N_enriched",     figs_dir = "figs/")

# Developmental time plots (requires a 'time' column in markers)
plot_mean_time( markerResList, thresh = 2, label = "MeanTime",       figs_dir = "figs/")
plot_order_comp(markerResList, thresh = 2,                           figs_dir = "figs/")
```

## Function reference

| Function | Description |
|---|---|
| `marker_enrich_table(res, markers)` | Score all cell types for one condition |
| `validate_markers(markers)` | Check markers input format |
| `split_markers(df, col, sep)` | Convert comma-separated string column to list-column |
| `plot_layer(res, name, thresh, filter)` | Dot plot of mean LFC per cell type; `filter=FALSE` shows all, `filter=TRUE` shows above-threshold only |
| `plot_comp(ordered_res_list, thresh, filter)` | Cross-condition dot plot; `filter=TRUE` anchors on above-threshold cell types, `filter=FALSE` shows all |
| `plot_comp_multi(res_list, thresh, filter)` | Run `plot_comp` with each condition as anchor |
| `plot_uniqueness(ordered_res_list, top)` | Top N cell types by uniqueness score (anchor mean − max of others) |
| `plot_uniqueness_multi(res_list, top)` | Run `plot_uniqueness` with each condition as anchor |
| `compute_auc(res_list, thresh)` | Return named vector of AUC scores |
| `plot_auc_curve(res_list, thresh, label)` | Smooth rank curves with AUC in caption |
| `plot_auc_bar(res_list, thresh, label)` | Bar plot of AUC per condition |
| `plot_n_enriched(res_list, thresh, label)` | Bar plot of enriched cell type count per condition |
| `plot_mean_time(res_list, thresh, label)` | Bar chart of mean developmental time of enriched cell types per condition |
| `plot_order_comp(res_list, thresh)` | Scatter of newly enriched cell types at each consecutive condition step, by developmental time |

All plot functions accept `figs_dir = NULL` to return the plot object without saving.

## Comparison to related methods

CTEnrich scores cell types by **mean LFC of their marker genes** — a simpler statistic than GSEA's running-sum enrichment score or GSVA's kernel-based score, but one that is on the same scale across all conditions without normalization. This makes cross-condition comparisons and AUC summaries directly interpretable as fold-change magnitudes.

- **GSEA / fgsea**: captures rank position of markers relative to all genes; NES is normalized per run and not cross-comparable across conditions.
- **camera (limma)**: competitive gene set test with p-values; designed for single contrasts, not multi-condition visualization.
- **GSVA**: per-sample scores from expression matrices; not designed for LFC input or cross-condition comparison.
- **AUCell**: single-cell focused; not applicable to bulk LFC results.

## References

Qiu C, et al. (2024) A single-cell time-lapse of mouse development from gastrulation to birth. *Nature*.
