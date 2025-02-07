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

add_marker_genes.over_representation <- function(meta_clusters, data.iteration, params, pvalue_threshold=0.001) {
  #' Get marker genes characteristic of each meta-cluster by conducting an over-representation test of the expressed genes.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
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
  cells_in_dataset <- ncol(data.iteration$expression.init)
  occurrences_of_genes.dataset <- table(data.iteration$ranking_of_genes.init)

  get_marker_genes.over_representation.cluster <- function(cluster) {
    cells_in_cluster <- length(cluster$cells)
    occurrences_of_genes.cluster <- table(data.iteration$ranking_of_genes.init[, cluster$cells])

    test_over_representation.gene <- function(gene) {
      test_over_representation(occurrences_of_genes.cluster[gene], cells_in_cluster,
                               occurrences_of_genes.dataset[gene], cells_in_dataset)}

    pvalues <- sapply(X=names(occurrences_of_genes.cluster), FUN=test_over_representation.gene, USE.NAMES=FALSE)
    adjusted_pvalues <- stats::p.adjust(pvalues, method="BH")
    marker_genes <- pvalues[adjusted_pvalues < pvalue_threshold]
    marker_genes <- sort(marker_genes)
    marker_genes <- -log10(marker_genes)
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
  markers <- sapply(X=meta_clusters, FUN="[[", "markers")
  markers <- unlist(unname(markers))
  markers <- table(names(markers))
  general_markers <- markers[markers > 1]

  get_specific_markers.cluster <- function(cluster) {setdiff(names(cluster$markers), names(general_markers))}
  for (i in 1:length(meta_clusters)) {
    meta_clusters[[i]][["specific_markers"]] <- get_specific_markers.cluster(meta_clusters[[i]])}
  return(meta_clusters)
}

get_characterized_clusters.markers_threshold <- function(meta_clusters, params, markers_threshold=22) {
  #' Check which meta-clusters are characterized, i.e. which meta-clusters have a minimum number of marker genes.
  #'
  #' If every meta-cluster is characterized, return them. Else, return none.
  #'
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param markers_threshold the minimal number of marker genes expected in a characterized cluster.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  #' @export
  #'
  is_characterized <- function(cluster) {length(cluster$markers) > markers_threshold}
  characterized_clusters <- Filter(f=is_characterized, x=meta_clusters)
  if (length(characterized_clusters) == length(meta_clusters)) {return(meta_clusters)}
  return(list())
}

get_characterized_clusters <- function(population, meta_clusters, data.iteration, params, figures) {
  #' Get characterized meta-clusters. They correspond to the sub-populations predicted by the clustering iteration.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param meta_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #'
  #' @return a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  #' @import glue
  #' @import grDevices
  #'
  #' @export
  #'
  meta_clusters <- params$marker_genes_strategy(meta_clusters, data.iteration, params)
  meta_clusters <- get_specific_markers(meta_clusters)
  if (figures) {
    plot <- draw_marker_genes(meta_clusters)
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_marker_genes.pdf"))
    print(plot)
    grDevices::dev.off()}
  characterized_clusters <- params$characterized_clusters_strategy(meta_clusters, params)
  return(characterized_clusters)
}

generate_color_scale <- function(labels){
  #' Generate a vector of colors equal to the number of identities in the sample.
  #'
  #' This function is directly copied from the repository of the SCpubr package.
  #' c.f. https://github.com/enblacar/SCpubr/blob/main/R/utils.R
  #'
  #' @param labels a vector of cluster labels.
  #'
  #' @return a named vector of colors.
  #'
  #' @import colorspace
  #' @import grDevices
  #'
  colors <- colorspace::qualitative_hcl(length(labels), palette = "Dark 3")
  colors <- grDevices::col2rgb(colors)
  colors <- grDevices::rgb2hsv(colors)
  colors["v", ] <- colors["v", ] - 0.1
  colors["s", ] <- colors["s", ] + 0.2
  colors["s", ][colors["s", ] > 1] <- 1
  colors <- grDevices::hsv(h = colors["h", ],
                           s = colors["s", ],
                           v = colors["v", ],
                           alpha = 1)
  names(colors) <- labels
  return(colors)
}

draw_marker_genes <- function(meta_clusters) {
  #' Get an upset-plot representing the marker genes predicted in each meta-cluster.
  #'
  #' @param meta_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  #' @return a plot.
  #'
  #' @import ggVennDiagram
  #' @import SCpubr
  #'
  #' @export
  #'
  get_markers.cluster <- function(cluster) {names(cluster$markers)}
  markers <- sapply(X=meta_clusters, FUN=get_markers.cluster)
  labels <- sapply(X=meta_clusters, FUN="[[", "label")
  colormap <- generate_color_scale(labels)

  # these functions are used to improve the aesthetics of the plots_____________
  add_grid_colors <- function(plot) {
    grid_colors <- rep(NA, length(plot[[1]]$data$name))
    grid_colors[1:length(labels)] <- labels
    plot[[1]]$layers[[1]]$aes_params$colour <- NULL
    plot[[1]]$layers[[2]]$aes_params$colour <- NULL
    plot[[1]] <- plot[[1]] + aes(colour=grid_colors) + scale_colour_manual(values=colormap)
    plot[[1]] <- plot[[1]] + theme(legend.position="none")
    return(plot)}

  add_intersection_colors <- function(plot) {
    intersection_colors <- rep(NA, length(plot[[2]]$data$name))
    intersection_colors[1:length(labels)] <- labels
    plot[[2]] <- plot[[2]] + aes(fill=intersection_colors) + scale_fill_manual(values=colormap)
    plot[[2]] <- plot[[2]] + scale_y_continuous(expand=expansion(mult=c(0, .05)))
    plot[[2]] <- plot[[2]] + theme(axis.text.y=element_text(vjust=0.25), legend.position="none")
    return(plot)}

  add_size_colors <- function(plot) {
    plot[[3]] <- plot[[3]] + aes(fill=labels) + scale_fill_manual(values=colormap)
    plot[[3]] <- plot[[3]] + theme(legend.position="none")
    return(plot)}
  #_____________________________________________________________________________

  plot <- ggVennDiagram::ggVennDiagram(markers, category.names=labels, force_upset=TRUE,
                                       order.set.by="name", order.intersect.by="none",
                                       relative_height=1.6, relative_width=0.4)
  plot <- add_grid_colors(plot)
  plot <- add_intersection_colors(plot)
  plot <- add_size_colors(plot)
  return(plot)
}
