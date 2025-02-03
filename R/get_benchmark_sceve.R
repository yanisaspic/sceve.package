"Functions called to benchmark an instance of the scEVE framework.

	2025/02/03 @yanisaspic"

get_extrinsic_clustering_metrics <- function(ground_truth, labels) {
  #' Using extrinsic metrics, measure the clustering performance.
  #' Extrinsic metrics compare the clustering results of an algorithm 
  #' to the ground truth of the authors of the dataset.
  #' The extrinsic metrics used are the Adjusted Rand index (ARI) and the Normalized Mutual Information (NMI).
  #' 
  #' @param ground_truth a named factor associating cells to their cluster annotations.
  #' @param labels a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector of numerics with two names: `ARI` and `NMI`.
  #' 
  #' @import aricode
  #'
  ground_truth <- ground_truth[names(labels)]
  extrinsic_clustering_metrics <- c("ARI"=aricode::ARI(ground_truth, labels),
                                    "NMI"=aricode::NMI(ground_truth, labels))
  return(extrinsic_clustering_metrics)
}

get_intrinsic_clustering_metrics <- function(expression.init, labels) {
  #' Using intrinsic metrics, measure the clustering performance.
  #' Intrinsic metrics compare the clustering results of an algorithm
  #' to the similarity in expression of the cells in a cluster,
  #' with regards to the dissimilarity in expression of the cells in different clusters.
  #' The intrinsic metrics used are the average neighborhood purity and the average silhouette index of the clusters.
  #' 
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param labels a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector of numerics with two names: `neighborhood_purity` and `silhouette_index`.
  #' 
  #' @import bluster
  #' 
  data <- t(expression.init)
  data <- data[names(labels), ]
  #neighborhood_purity <- bluster::neighborPurity(data, labels)
  silhouette_index <- bluster::approxSilhouette(data, labels)
  intrinsic_clustering_metrics <- c(#"neighborhood_purity"=mean(neighborhood_purity$purity),
                                    "silhouette_index"=mean(silhouette_index$width))
  return(intrinsic_clustering_metrics)
}

get_benchmark_sceve <- function(data, params) {
  #' Using both extrinsic and intrisic metrics, measure the clustering performance.
  #' 
  #' @param data a named list with two elements: `expression.init` and `ground_truth`.
  #' `expression.init` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' 
  #' @return a named vector of numerics with six names: `time`, `peak_memory_usage`, 
  #' `ARI`, `NMI`, `neighbordhood_purity` and `silhouette_index`.
  #' 
  memory <- gc(reset=TRUE)
  memory_usage0 <- memory[[11]] + memory[[12]]
  time0 <- Sys.time() #_________________________________________________________
  results <- sceve(data$expression.init, params, figures=FALSE, sheets=FALSE)
  time <- Sys.time()  #_________________________________________________________
  memory <- gc()
  memory_usage <- memory[[11]] + memory[[12]]
  
  metrics <- c("time"=time-time0, "peak_memory_usage"=memory_usage-memory_usage0)
  return(metrics)
}