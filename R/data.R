# markers_human: cell type marker gene sets from Qiu et al. 2024,
# converted from mouse to human orthologs via Ensembl BioMart.
#
# A data frame with 252 rows and 4 columns:
#   cell_type    - character; unique cell type name
#   marker_genes - list-column of character vectors; HGNC human gene symbols
#   germ_layer   - character; one of Ectoderm, Endoderm, Extraembryonic,
#                  Mesoderm, Other, PGC
#   time         - numeric; mouse embryonic day (e.g. 8.5 = E8.5);
#                  NA for cell types without a time assignment
#
# Source: Qiu C et al. (2024) A single-cell time-lapse of mouse development
#   from gastrulation to birth. Nature.
#   Orthologs mapped via Ensembl BioMart (GRCh38 / GRCm39).
#   Time annotations from typetime_df.rds (see data-raw/make_markers.R).
#
# Usage:
#   data(markers_human)
#   tbl <- marker_enrich_table(res, markers_human)
"markers_human"

# markers_mouse: cell type marker gene sets from Qiu et al. 2024,
# using original MGI mouse gene symbols.
#
# A data frame with 252 rows and 4 columns:
#   cell_type    - character; unique cell type name
#   marker_genes - list-column of character vectors; MGI mouse gene symbols
#   germ_layer   - character; one of Ectoderm, Endoderm, Extraembryonic,
#                  Mesoderm, Other, PGC
#   time         - numeric; mouse embryonic day (e.g. 8.5 = E8.5);
#                  NA for cell types without a time assignment
#
# Source: Qiu C et al. (2024) A single-cell time-lapse of mouse development
#   from gastrulation to birth. Nature.
#   Time annotations from typetime_df.rds (see data-raw/make_markers.R).
#
# Usage:
#   data(markers_mouse)
#   tbl <- marker_enrich_table(res, markers_mouse)
"markers_mouse"
