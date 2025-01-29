"Function called to load a scRNA-seq dataset presented in the vignette.

	2025/01/16 @yanisaspic"

get_cell_ids <- function(labels) {
  #' Get a vector of unique cell ids according to the label of each cell.
  #'
  #' @param labels a vector of cell types.
  #'
  #' @return a vector of characters.
  #'
  #' @import glue
  #'
  cell_ids <- c()
  counter <- list()
  for (cell_type in unique(labels)) {counter[[cell_type]] <- 1}
  for (n in 1:length(labels)) {
    cell_type <- labels[[n]]
    cell_ids[n] <- glue::glue("{cell_type}_{counter[[cell_type]]}")
    counter[[cell_type]] <- counter[[cell_type]] + 1}
  return(cell_ids)
}

get_data <- function(dataset) {
  #' Get a data.frame and the ground truth of a dataset loaded with TMExplorer.
  #'
  #' @param dataset a list loaded with TMExplorer::queryTME.
  #'
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @import SingleCellExperiment
  #' @import stats
  #' @import SummarizedExperiment
  #'
  throwaway_labels <- c("?", "NA", "", "Doublets")
  predictions <- SummarizedExperiment::colData(dataset)
  ground_truth <- stats::setNames(predictions$label, rownames(predictions))
  ground_truth <- ground_truth[!ground_truth %in% throwaway_labels]
  expression.init <- SingleCellExperiment::counts(dataset)
  expression.init <- expression.init[, names(ground_truth)]
  data <- list(expression.init=expression.init, ground_truth=ground_truth)
  return(data)
}

load_data <- function() {
  #' Load the scRNA-seq dataset used for the vignette.
  #' It corresponds to tissues within and surrounding human glioblastoma tumors.
  #' c.f. 10.1016/j.celrep.2017.10.030
  #' Darmanis, Spyros, et al.
  #' "Single-cell RNA-seq analysis of infiltrating neoplastic cells at the migrating front
  #' of human glioblastoma."
  #' Cell reports 21.5 (2017): 1399-1410.
  #'
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @import TMExplorer
  #'
  #' @export
  #'
  dataset <- TMExplorer::queryTME("GSE84465")[[1]]
  data <- get_data(dataset)
  cell_ids <- get_cell_ids(data$ground_truth)
  names(data$ground_truth) <- cell_ids
  colnames(data$expression.init) <- cell_ids
  return(data)
}
