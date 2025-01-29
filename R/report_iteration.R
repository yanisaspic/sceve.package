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
  for (drawing in c("extract_data", "base_clusters", "meta_clusters", "marker_genes")) {
    filename <- glue("{params$figures_path}/{population}_{cat}.pdf")
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

get_cell_memberships.binary_membership <- function(characterized_clusters, data.iteration, params) {
  #' Get a table reporting the cell composition of the characterized clusters.
  #' By calling this function, the cells assigned to the leftover cluster are hard clustered,
  #' i.e. their membership to the leftover is 1 and their membership to the other clusters is 0.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `ranking_of_genes` and `expression.init`.
  #' The three first elements correspond to the scRNA-seq expression matrix of a specific cell population,
  #' as well as the SeuratObject and the ranking of genes generated from this matrix.
  #' The last element corresponds to the full scRNA-seq expression matrix.
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
  data <- utils::stack(cells) %>% dplyr::rename(cells=.data$values, clusters=.data$ind)
  cells.iteration <- table(data$cells, data$clusters)
  cells.iteration <- as.data.frame.matrix(cells.iteration)
  return(cells.iteration)
}

report_cells <- function(characterized_clusters, data.iteration, records, params) {
  #' Save the cell composition of the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `ranking_of_genes` and `expression.init`.
  #' The three first elements correspond to the scRNA-seq expression matrix of a specific cell population,
  #' as well as the SeuratObject and the ranking of genes generated from this matrix.
  #' The last element corresponds to the full scRNA-seq expression matrix.
  #' @param records a named list, with three data.frames: `cells`, `markers` and `meta`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are characterized clusters, and cell memberships are reported in the table.
  #'
  cell_memberships <- params$cell_memberships_strategy(characterized_clusters, data.iteration, params)
  records$cells <- cbind(records$cells, cell_memberships)
  return(records$cells)
}

report_metadata <- function(population, characterized_clusters, records) {
  #' Report the metadata regarding characterized clusters in the corresponding records sheet.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param records a named list, with three data.frames: `cells`, `markers` and `meta`.
  #'
  #' @return a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  get_metadata.cluster <- function(cluster) {
    cells_of_population <- get_cells_of_population(cluster$label, records$cells)
    c(size=length(cells_of_population), robustness=cluster$robustness,
      parent=population, clustering_status="PENDING")}
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
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `ranking_of_genes` and `expression.init`.
  #' The three first elements correspond to the scRNA-seq expression matrix of a specific cell population,
  #' as well as the SeuratObject and the ranking of genes generated from this matrix.
  #' The last element corresponds to the full scRNA-seq expression matrix.
  #' @param records a named list, with three data.frames: `cells`, `markers` and `meta`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named list, with three data.frames: `cells`, `markers` and `meta`.
  #'
  records$cells <- report_cells(characterized_clusters, data.iteration, records, params)
  records$meta <- report_metadata(population, characterized_clusters, records)
  return(records)
}
