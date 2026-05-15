# markers_human: cell type marker gene sets from Qiu et al. 2024,
# converted from mouse to human orthologs via Ensembl BioMart.
#
# A data frame with 252 rows and 3 columns:
#   cell_type    - character; unique cell type name
#   marker_genes - list-column of character vectors; HGNC human gene symbols
#   germ_layer   - character; one of Ectoderm, Endoderm, Extraembryonic,
#                  Mesoderm, Other, PGC
#
# Source: Qiu C et al. (2024) A single-cell time-lapse of mouse development
#   from gastrulation to birth. Nature.
#   Orthologs mapped via Ensembl BioMart (GRCh38 / GRCm39).
#
# Usage:
#   data(markers_human)
#   tbl <- marker_enrich_table(res, markers_human)
"markers_human"
