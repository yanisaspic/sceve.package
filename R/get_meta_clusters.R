"Functions used to get meta-clusters, i.e. robust clusters and a leftover cluster.

	2025/01/21 @yanisaspic"

get_transaction_database <- function(base_clusters) {
  #' Get a transaction database from the base clusters predicted.
  #'
  #' @param base_clusters a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @return a transactions object, where every cell is a transaction and every predicted
  #' cluster is an item.
  #'
  #' @import arules
  #' @import glue
  #' @import rlang
  #' @import utils
  #'
  tmp <- rlang::hash(base_clusters) # use hash to prevent issues with parallel calculations.
  path <- glue::glue("./{tmp}.tmp")
  utils::write.table(base_clusters, file=path, col.names=FALSE, row.names=FALSE)
  transaction_database <- arules::read.transactions(path)
  file.remove(path)
  return(transaction_database)
}

get_associations <- function(transaction_database) {
  #' Measure the strength of the association between pairs of predicted clusters,
  #' according to the confidence metric of the frequent itemset mining framework.
  #' The confidence(A->C) corresponds to the proportion of cells in the cluster A, that
  #' are also in the cluster C. Note that confidence(A->C) can differ from confidence(C->A).
  #' Associations involving a minority of cells (i.e. confidence < 0.5) are filtered out.
  #'
  #' @param transaction_database a transactions object, where every cell is a transaction and every predicted
  #' cluster is an item.
  #'
  #' @return a data.frame associating clusters `A` and `C` to their `confidence(A->C)`.
  #'
  #' @import arules
  #'
  data <- arules::apriori(transaction_database,
                          parameter = list(support=0.001, confidence=0.5, minlen=2, maxlen=2, target="rules"))
  associations <- arules::DATAFRAME(data)

  get_item <- function(x) {substr(x, 2, nchar(x)-1)}
  associations <- list(A=sapply(X=as.character(associations$LHS), FUN=get_item),
                                C=sapply(X=as.character(associations$RHS), FUN=get_item),
                                confidence=as.numeric(associations$confidence))
  associations <- do.call(cbind, associations)
  associations <- data.frame(associations, row.names = NULL)
  return(associations)
}

get_strong_similarities <- function(associations) {
  #' Get pairs of clusters that share a strong similarity.
  #' A strong similarity is observed when two clusters share the majority of their cells,
  #' i.e. confidence(A->C) > 0.5 & confidence(C->A) > 0.5.
  #' It corresponds to the minimum proportion of cells shared between two clusters.
  #'
  #' @param associations a data.frame associating two clusters `A` and `C`, to their `confidence(A->C)`.
  #'
  #' @return a data.frame associating two clusters `A` and `C`, to their `similarity`.
  #'
  associations <- associations[order(associations$confidence, decreasing = TRUE),]

  get_id.row <- function(row) {
    elements <- sort(row[c("A", "C")])
    id <- paste(elements, collapse=".")}
  ids <- apply(X=associations, MARGIN=1, FUN=get_id.row)

  strong_similarities <- associations[duplicated(ids),]
  colnames(strong_similarities)[3] <- "similarity"
  return(strong_similarities)
}

get_cells_of_subgraph <- function(subgraph, base_clusters) {
  #' Get the cells at the intersection of all the clusters in a subgraph.
  #'
  #' @param subgraph a subgraph corresponding to strongly similar clusters.
  #' @param base_clusters a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @return a vector of cells.
  #'
  clusters_of_subgraph <- names(igraph::V(subgraph))
  cluster_is_in_subgraph <- function(cluster) {cluster %in% clusters_of_subgraph}
  data <- apply(X=base_clusters, MARGIN=c(1,2), FUN=cluster_is_in_subgraph)
  cells <- rownames(data[rowSums(data)==length(clusters_of_subgraph), ])
  return(cells)
}

get_subgraphs <- function(population, strong_similarities, base_clusters, params) {
  #' Model a graph where every node is a base cluster predicted, and edges are strong similarities.
  #' The graph is disjoint, and subgraphs (i.e. connected components) are extracted from it.
  #' They will be used to identify robust clusters.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param strong_similarities a data.frame associating two clusters `A` and `C`, to their `similarity`.
  #' @param base_clusters a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a list where every element is a pool of cells grouped together by multiple clustering methods.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @import igraph
  #' @import glue
  #'
  data <- igraph::graph.data.frame(strong_similarities)
  subgraphs <- igraph::decompose.graph(data)
  n_methods <- ncol(base_clusters)
  theoretical_robustness <- n_methods * (n_methods - 1) / 2

  # these functions are used to format the subgraph data in a list______________
  get_clusters_of_subgraph <- function(subgraph) {names(igraph::V(subgraph))}
  get_robustness_of_subgraph <- function(subgraph) {
    sum(as.numeric(igraph::E(subgraph)$similarity)) / theoretical_robustness}
  get_method_of_renamed_cluster <- function(renamed_cluster) {
    strsplit(renamed_cluster, split="_")[[1]][1]}
  get_methods_of_subgraph <- function(subgraph) {
    unique(sapply(X=get_clusters_of_subgraph(subgraph), FUN=get_method_of_renamed_cluster))}

  get_subgraph <- function(subgraph) {
    list(base_clusters=get_clusters_of_subgraph(subgraph),
         robustness=get_robustness_of_subgraph(subgraph),
         clustering_methods=get_methods_of_subgraph(subgraph),
         cells=get_cells_of_subgraph(subgraph, base_clusters))}
  #_____________________________________________________________________________

  subgraphs <- lapply(X=subgraphs, FUN=get_subgraph)
  subgraph_has_intersection <- function(subgraph) {length(subgraph$cells) > 0}
  subgraphs <- Filter(f=subgraph_has_intersection, x=subgraphs)
  robustnesses <- sapply(subgraphs, "[[", "robustness")
  subgraphs <- subgraphs[order(-robustnesses)]
  for (i in 1:length(subgraphs)) {subgraphs[[i]][["label"]] <- glue::glue("{population}.{i}")}
  return(subgraphs)
}

add_leftover_cluster <- function(population, robust_clusters, data.iteration) {
  #' Add a leftover cluster to the list of robust clusters.
  #' It corresponds to a group of cells unassigned to any robust cluster.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param robust_clusters a list where every element is a pool of cells grouped together by multiple clustering methods.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @import glue
  #'
  n <- length(robust_clusters) + 1
  cells_of_population <- colnames(data.iteration$expression)
  cells_in_robust_clusters <- unlist(sapply(robust_clusters, "[[", "cells"))
  leftover_cells <- setdiff(cells_of_population, cells_in_robust_clusters)
  if (length(leftover_cells) > 0) {
    leftover_cluster <- list(base_clusters=c(), clustering_methods=c(), robustness=0,
                             cells=leftover_cells, label=glue::glue("{population}.{n}"))
    robust_clusters[[n]] <- leftover_cluster}
  return(robust_clusters)
}

get_meta_clusters <- function(population, base_clusters, data.iteration, records, params, figures) {
  #' Extract robust clusters and a leftover cluster from a set of base clusters predicted
  #' with multiple clustering methods.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param base_clusters a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @import glue
  #' @import grDevices
  #' @import stats
  #'
  #' @export
  #'
  transaction_database <- get_transaction_database(base_clusters)
  associations <- get_associations(transaction_database)
  strong_similarities <- get_strong_similarities(associations)
  subgraphs <- get_subgraphs(population, strong_similarities, base_clusters, params)
  robustness_threshold <- max(params$robustness_threshold, records$meta[population, "robustness"])
  subgraph_is_robust_cluster <- function(subgraph) {subgraph$robustness > params$robustness_threshold}
  robust_clusters <- Filter(f=subgraph_is_robust_cluster, x=subgraphs)
  meta_clusters <- add_leftover_cluster(population, robust_clusters, data.iteration)
  meta_clusters <- stats::setNames(meta_clusters, sapply(X=meta_clusters, FUN="[[", "label"))

  if (figures) {
    plot <- draw_meta_clusters(meta_clusters, data.iteration)
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_meta_clusters.pdf"))
    print(plot)
    grDevices::dev.off()}
  return(meta_clusters)
}

draw_meta_clusters <- function(meta_clusters, data.iteration) {
  #' Get composite U-MAP plots representing the meta-clusters predicted with the scEVE framework.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #'
  #' @return a plot.
  #'
  #' @import ggplot2
  #' @import glue
  #' @import gridExtra
  #' @import SCpubr
  #' @import stats
  #'
  #' @export
  #'
  get_plot_coloring.cluster <- function(cluster) {
    plot_label <- glue::glue("{cluster$label} ({round(cluster$robustness, 2)})")
    stats::setNames(rep(plot_label, length(cluster$cells)), cluster$cells)}
  plot_coloring <- sapply(X=meta_clusters, FUN=get_plot_coloring.cluster)
  plot_coloring <- unlist(unname(plot_coloring))
  data.iteration$SeuratObject[["meta_clusters"]] <- as.factor(plot_coloring)

  plot <- SCpubr::do_DimPlot(data.iteration$SeuratObject, split.by="meta_clusters", legend.position="none")
  for (i in 1:length(plot)) {
    plot[[i]][[1]] <- plot[[i]][[1]] +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title=element_text(hjust=0.5, margin=margin(1, 0, 0, 0)),
                     panel.background=element_rect(fill="lightgrey"),
                     legend.position="none", axis.title=element_blank())}
  return(plot)
}
