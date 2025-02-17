"Functions called to benchmark an instance of the scEVE framework.

	2025/02/03 @yanisaspic"

get_data.bluster <- function(expression.init) {
  #' Get a matrix usable for the functions of bluster.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return an object associating selected genes and cells to their reduced dimensions.
  #'
  #' @import scater
  #' @import scran
  #' @import scuttle
  #' @import SingleCellExperiment
  #'
  data <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(expression.init)))
  data <- scuttle::logNormCounts(data)
  variable_genes <- scran::getTopHVGs(scran::modelGeneVar(data), n=5000)
  set.seed(1)
  data <- scater::runPCA(data, ncomponents=20, subset_row=variable_genes)
  data <- SingleCellExperiment::reducedDim(data)
  return(data)
}

get_clustering_metrics.intrinsic <- function(expression.init, preds) {
  #' Using intrinsic metrics, measure a clustering performance.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param preds a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector with one name: `SI`.
  #'
  #' @import bluster
  #'
  n_clusters_predicted <- length(unique(preds))
  if (n_clusters_predicted < 2) {return(c("SI"=NA))}
  data <- get_data.bluster(expression.init)
  silhouette_index <- bluster::approxSilhouette(data, preds)
  clustering_metrics.intrinsic <- c("SI"=mean(silhouette_index$width))
  return(clustering_metrics.intrinsic)
}

get_clustering_metrics <- function(data, preds) {
  #' Using both intrinsic and extrinsic clustering metrics, measure a clustering performance.
  #'
  #' Extrinsic metrics compare cluster predictions to the cell annotations of the dataset.
  #' Intrinsic metrics compare the gene expression of cells in and out of their clusters.
  #' The `ARI`, the `NMI` and the `Purity` are extrinsic metrics.
  #' The `SI` is an intrinsic metric.
  #' For every metric, higher is better and the maximum value is 1.
  #'
  #' @param data a named list with two elements: `expression.init` and `ground_truth`.
  #' `expression.init` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param preds a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector with four names: `ARI`, `NMI`, `Purity` and `SI`.
  #'
  #' @import aricode
  #' @import funtimes
  #'
  if (length(preds) < 2) {return(c("ARI"=NA, "NMI"=NA, "Purity"=NA, "SI"=NA))}
  ground_truth <- data$ground_truth[names(preds)]
  expression.init <- data$expression.init[, names(preds)]
  clustering_metrics <- c("ARI"=aricode::ARI(ground_truth, preds),
                          "NMI"=aricode::NMI(ground_truth, preds),
                          "Purity"=funtimes::purity(ground_truth, preds)$pur,
                          get_clustering_metrics.intrinsic(expression.init, preds))
  return(clustering_metrics)
}

get_benchmark_sceve.data <- function(data, params, method_label) {
  #' Using computational as well as intrinsic and extrinsic clustering metrics, measure
  #' the performance of an instance of the scEVE framework on a dataset.
  #'
  #' @param data a named list with two elements: `expression.init` and `ground_truth`.
  #' `expression.init` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param method_label a character.
  #'
  #' @return a data.frame with eight columns: `method`, `time (s)`,
  #' `peak_memory_usage (Mb)`, `ARI`, `NMI`, `Purity`, `SI` and `n_cells`.
  #'
  #' @import glue
  #'
  #' @export
  #'
  get_memory_usage <- function(memory) {memory[[11]] + memory[[12]]}
  memory_usage.init <- get_memory_usage(gc(reset=TRUE))
  time.init <- Sys.time()
  results <- sceve(data$expression.init, params, figures=FALSE, sheets=FALSE)
  time <- as.numeric(Sys.time() - time.init, units="secs")
  peak_memory_usage <- get_memory_usage(gc()) - memory_usage.init

  population_is_leftover <- function(population) {
    (results$records$meta[population, "robustness"] == 0)}
  leftover_populations <- Filter(population_is_leftover, levels(results$preds))

  benchmark <- c("method"=method_label, "time (s)"=time,
                 "peak_memory_usage (Mb)"=peak_memory_usage,
                 get_clustering_metrics(data, results$preds),
                 "n_cells"=length(results$preds))
  benchmark.star <- c("method"=glue::glue("{method_label}*"), "time (s)"=time,
                      "peak_memory_usage (Mb)"=peak_memory_usage,
                      get_clustering_metrics(data, results$preds[!results$preds %in% leftover_populations]),
                      "n_cells"=length(results$preds[!results$preds %in% leftover_populations]))
  ground_truth <- c("method"="ground_truth", "time (s)"=NA, "peak_memory_usage (Mb)"=NA,
                    get_clustering_metrics(data, data$ground_truth),
                    "n_cells"=length(data$ground_truth))
  benchmark_sceve <- as.data.frame(rbind(benchmark, benchmark.star, ground_truth))
  return(benchmark_sceve)
}

get_benchmark_sceve <- function(datasets, params, method_label) {
  #' Using computational as well as intrinsic and extrinsic clustering metrics, measure
  #' the performance of an instance of the scEVE framework on multiple datasets.
  #'
  #' @param datasets a vector of datasets (cf. `sceve::get_datasets()`).
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param method_label a character.
  #'
  #' @return a data.frame with nine columns: `method`, `dataset`, `time (s)`,
  #' `peak_memory_usage (Mb)`, `ARI`, `NMI`, `Purity`, `SI` and `n_cells`.
  #'
  #' @export
  #'
  get_benchmark_sceve.dataset <- function(dataset) {
    data <- load_data(dataset)
    benchmark <- get_benchmark_sceve.data(data, params, method_label)
    benchmark[, "dataset"] <- dataset
    return(benchmark)}

  benchmarks <- lapply(X=datasets, FUN=get_benchmark_sceve.dataset)
  benchmark <- do.call(rbind, benchmarks)
  rownames(benchmark) <- NULL
  return(benchmark)
}
