# CTEnrich

Cell type enrichment scoring and visualization for bulk RNA-seq.

Each condition's DESeq2 result is scored against a set of cell-type marker gene lists by averaging the LFC of each cell type's markers. High mean LFC = that cell type's signature is upregulated relative to the baseline.

> **Recommended baseline:** pluripotent cells (iPSCs or ESCs). Same-experiment, same-library-prep samples are best; public transcriptomes may introduce batch effects.

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
| `cell_type` | yes | Unique cell type name |
| `marker_genes` | yes | List-column of character vectors (gene symbols) |
| `germ_layer` | no | Used for colour in plots |
| `time` | no | Developmental time; required for `plot_mean_time` and `plot_order_comp` |

Marker sets should come from scRNA-seq data matching your biological context (e.g. top DE genes per cluster). Gene symbols must match `rownames(res)` (both HGNC or both MGI).

Use `split_markers()` to convert a comma-separated string column to a list-column, and `validate_markers()` to check the input. Extra columns pass through `marker_enrich_table()` unchanged.

## Built-in marker sets

Two marker sets from the Qiu et al. (2024) mouse developmental atlas: 252 cell types, six germ layers, ~19 markers each.

| Dataset | Symbols | Use with |
|---|---|---|
| `markers_human` | HGNC (mouse orthologs via Ensembl BioMart) | Human bulk RNA-seq |
| `markers_mouse` | MGI (original mouse symbols) | Mouse bulk RNA-seq |

```r
data(markers_human)
data(markers_mouse)
```

## Quick start

```r
library(CTEnrich)

data(markers_human)

# Or supply your own
markers <- data.frame(
  cell_type    = c("Cardiomyocyte", "Hepatocyte", "Neuron"),
  marker_genes = c("TNNT2, MYH7, ACTC1", "ALB, AFP, HNF4A", "MAP2, SYP, NCAM1"),
  germ_layer   = c("Mesoderm", "Endoderm", "Ectoderm")
)
markers <- split_markers(markers, "marker_genes", sep = ", ")
validate_markers(markers)

# Score (resList = named list of DESeq2 results)
markerResList <- lapply(resList, marker_enrich_table, markers = markers_human)

# Per-condition plots
plot_layer(markerResList[["condition_a"]], name = "condition_a", filter = TRUE)

# Cross-condition (filter = TRUE: above threshold; FALSE: all cell types)
plot_comp_multi(markerResList, thresh = 2, filter = TRUE,  figs_dir = "figs/")
plot_comp_multi(markerResList, thresh = 2, filter = FALSE, figs_dir = "figs/")

# Uniqueness scores
plot_uniqueness_multi(markerResList, top = 30, figs_dir = "figs/")

# Summary plots
plot_auc_curve( markerResList, thresh = 2, label = "AUC_comparison", figs_dir = "figs/")
plot_auc_bar(   markerResList, thresh = 2, label = "AUC_bar",        figs_dir = "figs/")
plot_n_enriched(markerResList, thresh = 2, label = "N_enriched",     figs_dir = "figs/")

# Developmental time plots (requires 'time' column in markers)
plot_mean_time( markerResList, thresh = 2, label = "MeanTime",       figs_dir = "figs/")
plot_order_comp(markerResList, thresh = 2,                           figs_dir = "figs/")
```

## Function reference

| Function | Description |
|---|---|
| `marker_enrich_table(res, markers)` | Score all cell types for one condition |
| `validate_markers(markers)` | Check markers input format |
| `split_markers(df, col, sep)` | Convert comma-separated string column to list-column |
| `plot_layer(res, name, thresh, filter)` | Mean LFC dot plot; `filter=FALSE` all cell types, `filter=TRUE` above threshold |
| `plot_comp(ordered_res_list, thresh, filter)` | Cross-condition dot plot; `filter=TRUE` anchors on above-threshold cell types |
| `plot_comp_multi(res_list, thresh, filter)` | Run `plot_comp` with each condition as anchor |
| `plot_uniqueness(ordered_res_list, top)` | Top N cell types by uniqueness score (anchor mean − max of others) |
| `plot_uniqueness_multi(res_list, top)` | Run `plot_uniqueness` with each condition as anchor |
| `compute_auc(res_list, thresh)` | Named vector of AUC scores |
| `plot_auc_curve(res_list, thresh, label)` | Smooth rank curves with AUC in caption |
| `plot_auc_bar(res_list, thresh, label)` | Bar plot of AUC per condition |
| `plot_n_enriched(res_list, thresh, label)` | Bar plot of enriched cell type count per condition |
| `plot_mean_time(res_list, thresh, label)` | Mean developmental time of enriched cell types per condition |
| `plot_order_comp(res_list, thresh)` | Newly enriched cell types at each consecutive condition step, by developmental time |

All plot functions accept `figs_dir = NULL` to return the plot object without saving.

## Comparison to related methods

CTEnrich scores by **mean LFC of marker genes**, simpler than GSEA's running-sum or GSVA's kernel score, but on the same scale across all conditions without normalization, making cross-condition comparisons interpretable directly as fold-change.

- **GSEA / fgsea**: NES is normalized per run, not cross-comparable across conditions.
- **camera (limma)**: single-contrast p-values, not built for multi-condition visualization.
- **GSVA**: per-sample scores from expression matrices; not designed for LFC input.
- **AUCell**: single-cell focused; not applicable to bulk LFC.

## References

Qiu C, et al. (2024) A single-cell time-lapse of mouse development from gastrulation to birth. *Nature*.
