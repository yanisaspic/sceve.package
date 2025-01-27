"Functions called to identify characterized clusters, i.e. clusters with sufficient marker genes.

	2025/01/24 @yanisaspic"

test_over_representation <- function(cells_in_cluster.expressing, cells_in_pool.expressing,
                                     cells_in_cluster, cells_in_pool) {
  #' Test the over-representation of a gene in a cluster with regards to the frequency of its expression
  #' in the cluster and in the complete pool of cells.
  #'
  #' @param cells_in_cluster.expressing the number of cells expressing the gene in the cluster.
  #' @param cells_in_pool.expressing the number of cells expressing the gene in the pool of cells.
  #' @param cells_in_cluster the number of cells in the cluster.
  #' @param cells_in_pool the number of cells in the pool of cells.
  #'
  #' @return the p-value of the over-representation test.
  #'
  #' @import stats
  #'
  q <- cells_in_cluster.expressing - 1 # phyper measures P(X > q), but we want P(X >= cells_in_cluster.expressing)
  n <- cells_in_pool - cells_in_pool.expressing
  pvalue <- stats::phyper(q, cells_in_pool.expressing, n, cells_in_cluster, lower.tail=FALSE)
  return(pvalue)
}

get_marker_genes.over_representation <- function(meta_clusters, data.iteration, params, pvalue_threshold=0.05) {
  #' Get marker genes characteristic of each meta-cluster by conducting an over-representation test of the expressed genes.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with three names: `expression`, `SeuratObject` and `ranking_of_genes`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population,
  #' as well as the SeuratObject and the ranking of genes generated from this matrix.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param pvalue_threshold the pvalue threshold below which genes are considered as marker genes.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' @import stats
  #'
  #' @export
  #'
  cells_in_iteration <- ncol(data.iteration$ranking_of_genes)
  occurrences_of_genes.iteration <- table(data.iteration$ranking_of_genes)

  get_marker_genes.over_representation.cluster <- function(cluster) {
    cells_in_cluster <- length(cluster$cells)
    occurrences_of_genes.cluster <- table(data.iteration$ranking_of_genes[, cluster$cells])

    test_over_representation.gene <- function(gene) {
      test_over_representation(occurrences_of_genes.cluster[gene], cells_in_cluster,
                               occurrences_of_genes.iteration[gene], cells_in_iteration)}

    pvalues <- sapply(X=names(occurrences_of_genes.cluster), FUN=test_over_representation.gene, USE.NAMES=FALSE)
    adjusted_pvalues <- stats::p.adjust(pvalues, method="BH")
    marker_genes <- pvalues[adjusted_pvalues < pvalue_threshold]
    marker_genes <- sort(marker_genes)
    return(marker_genes)
  }

  for (i in 1:length(meta_clusters)) {
    meta_clusters[[i]][["markers"]] <- get_marker_genes.over_representation.cluster(meta_clusters[[i]])}
  return(meta_clusters)
}

get_specific_markers <- function(meta_clusters) {
  #' Get marker genes specific to each meta-cluster.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  markers <- unlist(sapply(X=meta_clusters, FUN="[[", "markers"))
  markers <- table(names(markers))
  general_markers <- markers[markers > 1]

  get_specific_markers.cluster <- function(cluster) {setdiff(names(cluster$markers), names(general_markers))}
  for (i in 1:length(meta_clusters)) {
    meta_clusters[[i]][["specific_markers"]] <- get_specific_markers.cluster(meta_clusters[[i]])}
  return(meta_clusters)
}

get_characterized_clusters.specific_markers_threshold <- function(robust_clusters, params,
                                                                  specific_markers_threshold=22) {
  #' Get characterized robust clusters, i.e. robust clusters with a minimal number of specific marker genes.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param specific_markers_threshold the minimal number of specific marker genes expected in a characterized cluster.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  #' @export
  #'
  is_characterized <- function(cluster) {length(cluster$specific_markers) > specific_markers_threshold}
  characterized_clusters <- Filter(f=is_characterized, x=robust_clusters)
  return(characterized_clusters)
}
