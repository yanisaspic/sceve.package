"Functions called to report the results of a scEVE clustering iteration.

	2025/01/29 @yanisaspic"

get_drawings_paths <- function(population, params) {
  #' Get the paths leading to each informative figure drawn during the clustering iteration.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a vector of paths.
  #'
  drawings_paths <- c()
  for (drawing in c("extract_data", "base_clusters", "robust_clusters", "marker_genes")) {
    filename <- glue("{params$figures_path}/{population}_{drawing}.pdf")
    if (file.exists(filename)) {drawings_paths <- c(drawings_paths, filename)}}
  return(drawings_paths)
}

merge_drawings <- function(population, params) {
  #' Merge the informative figures drawn during the clustering iteration in a single .pdf file.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @import glue
  #' @import qpdf
  #'
  drawings_paths <- get_drawings_paths(population, params)
  qpdf::pdf_combine(input=drawings_paths, output=glue("{params$figures_path}/{population}.pdf"))
  unlink(drawings_paths)
}

get_cluster_memberships.binary_membership <- function(characterized_clusters, data.iteration, params) {
  #' Get a table reporting the cell composition of the characterized clusters.
  #'
  #' By calling this function, the cells assigned to the leftover cluster are hard clustered,
  #' i.e. their membership to the leftover is 1 and their membership to the other clusters is 0.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are characterized clusters, and cell memberships are reported in the table.
  #'
  #' @import dplyr
  #' @importFrom rlang .data
  #' @import utils
  #'
  #' @export
  #'
  cells <- sapply(X=characterized_clusters, FUN="[[", "cells")
  data <- utils::stack(cells) %>% dplyr::rename(cells=values, clusters=ind)
  cells.iteration <- table(data$cells, data$clusters)
  cells.iteration <- as.data.frame.matrix(cells.iteration)
  return(cells.iteration)
}

report_cells <- function(characterized_clusters, data.iteration, records, params) {
  #' Report the cell composition of the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are characterized clusters, and cell memberships are reported in the table.
  #'
  cluster_memberships <- params$cluster_memberships_strategy(characterized_clusters, data.iteration, params)
  records$cells <- merge(records$cells, cluster_memberships, by="row.names", all=TRUE)
  records$cells <- transform(records$cells, row.names=Row.names, Row.names=NULL)
  records$cells[is.na(records$cells)] <- 0
  return(records$cells)
}

report_markers <- function(characterized_clusters, records) {
  #' Report the marker genes of the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to their marker genes.
  #' Its rows are genes, its columns are predicted populations,
  #' and the strength of the characterization (e.g. log10-transformed pvalues) are reported in the table.
  #'
  #' @import stats
  #'
  get_column <- function(cluster) {
    column <- stats::setNames(data.frame(cluster$markers), cluster$label)
    column[setdiff(rownames(records$markers), rownames(column)),] <- 0
    return(column)}
  columns <- lapply(X=characterized_clusters, FUN=get_column)
  for (c in columns) {
    records$markers <- merge(records$markers, c, by="row.names", all=TRUE)
    records$markers <- transform(records$markers, row.names=Row.names, Row.names=NULL)}
  return(records$markers)
}

report_methods <-function(characterized_clusters, records) {
  #' Report the clustering methods used to predict the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to the clustering methods leveraged to predict them.
  #' Its rows are clustering methods, its columns are predicted populations, and binary values are reported in the table.
  #'
  methods <- sapply(X=characterized_clusters, FUN="[[", "clustering_methods")
  data <- utils::stack(methods) %>% dplyr::rename(methods=values, clusters=ind)
  methods.iteration <- table(data$methods, data$clusters)
  methods.iteration <- as.data.frame.matrix(methods.iteration)
  records$methods <- merge(records$methods, methods.iteration, by="row.names", all=TRUE)
  records$methods <- transform(records$methods, row.names=Row.names, Row.names=NULL)
  records$methods[is.na(records$methods)] <- 0
  return(records$methods)
}

report_metadata <- function(population, characterized_clusters, records) {
  #' Report the metadata regarding characterized clusters in the corresponding records sheet.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  get_metadata.cluster <- function(cluster) {
    cells_of_population <- get_cells_of_population(cluster$label, records$cells)
    metadata.cluster <- c(size=length(cells_of_population), robustness=cluster$robustness,
                          parent=population, clustering_status="PENDING")
    return(metadata.cluster)}
  metadata.iteration <- lapply(X=characterized_clusters, FUN=get_metadata.cluster)
  metadata.iteration <- do.call(rbind, metadata.iteration)
  records$meta <- rbind(records$meta, metadata.iteration)
  for (col in c("robustness", "size")) {records$meta[, col] <- as.numeric(records$meta[, col])}
  return(records$meta)
}

report_iteration <- function(population, characterized_clusters, data.iteration, records, params) {
  #' Report multiple information regarding the characterized clusters predicted during a clustering iteration.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param records a named list, with three data.frames: `cells`, `markers` and `meta`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @export
  #'
  records$cells <- report_cells(characterized_clusters, data.iteration, records, params)
  records$markers <- report_markers(characterized_clusters, records)
  records$meta <- report_metadata(population, characterized_clusters, records)
  records$methods <- report_methods(characterized_clusters, records)
  return(records)
}
