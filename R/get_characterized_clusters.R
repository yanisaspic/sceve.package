"Functions called to identify characterized clusters, i.e. clusters with sufficient marker genes.

	2025/01/24 @yanisaspic"

get_marker_genes.seurat_findmarkers <- function(cluster, data.iteration, params, FC_threshold=4, pvalue_threshold=0.001) {
  #' Get marker genes by calling the function `FindMarkers` of the Seurat package.
  #'
  #' A named vector associating marker genes to their log2FC is returned.
  #' Only genes with log2FC > 4 and p-values (corrected with Bonferonni) < 0.001 are returned.
  #'
  #' @param cluster a named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration a named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param FC_threshold a numeric.
  #' @param pvalue_threshold a numeric.
  #'
  #' @return a named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`,`markers` and `robustness`.
  #'
  #' @import Seurat
  #' @import stats
  #'
  cells_of_iteration <- colnames(data.iteration$expression)
  if (length(cells_of_iteration) == length(cluster$cells)) {return(c())}
  is_in_cluster <- function(cell) {ifelse(cell %in% cluster$cells, 1, 0)}
  Seurat::Idents(object=data.iteration$SeuratObject) <- factor(is_in_cluster(cells_of_iteration))
  markers <- Seurat::FindMarkers(data.iteration$SeuratObject, ident.1=1)
  markers <- markers[(markers$avg_log2FC > FC_threshold) & (markers$p_val_adj < pvalue_threshold), ]
  markers <- stats::setNames(markers$avg_log2FC, rownames(markers))
  markers <- markers[order(-markers)]
  return(markers)
}

add_leftover_cluster <- function(population, robust_clusters, data.iteration, params) {
  #' Add a leftover cluster to the list of robust clusters.
  #'
  #' It corresponds to a group of cells unassigned to any robust cluster.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param robust_clusters list where every element is a robust pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #' @param data.iteration a named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, and `robustness`.
  #'
  #' @import glue
  #' @import stats
  #'
  n_clusters <- length(robust_clusters)
  cells_of_population <- colnames(data.iteration$expression)
  cells_in_robust_clusters <- unlist(sapply(robust_clusters, "[[", "cells"))
  leftover_cells <- setdiff(cells_of_population, cells_in_robust_clusters)
  if (length(leftover_cells) > 0) {
    leftover_cluster <- list(base_clusters=c(), clustering_methods=c(), robustness=0,
                             cells=leftover_cells, label=glue::glue("{population}.L"))
    leftover_cluster[["markers"]] <- params$marker_genes_strategy(leftover_cluster, data.iteration, params)
    robust_clusters[[n_clusters + 1]] <- leftover_cluster}
  robust_clusters <- stats::setNames(robust_clusters, sapply(X=robust_clusters, FUN="[[", "label"))
  return(robust_clusters)
}

get_characterized_clusters.markers_threshold <- function(clusters, data.iteration, params, markers_threshold=10) {
  #' Filter out uncharacterized clusters, i.e. clusters with too little marker genes.
  #'
  #' @param clusters a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #' @param data.iteration a named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param markers_threshold the minimal number of marker genes expected in a characterized cluster.
  #'
  #' @return a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' @export
  #'
  is_characterized <- function(cluster) {length(cluster$markers) >= markers_threshold}
  characterized_clusters <- Filter(f=is_characterized, x=clusters)
  return(characterized_clusters)
}

get_characterized_clusters <- function(population, robust_clusters, data.iteration, params, figures) {
  #' Get characterized meta-clusters. They correspond to the sub-populations predicted by the clustering iteration.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param robust_clusters a list where every element is a pool of cells.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #' @param data.iteration named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #'
  #' @return a list where every element is a characterized pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' @import glue
  #' @import grDevices
  #'
  #' @export
  #'
  add_markers <- function(cluster) {
    cluster[["markers"]] <- params$marker_genes_strategy(cluster, data.iteration, params)
    return(cluster)}
  robust_clusters <- lapply(X=robust_clusters, FUN=add_markers)
  characterized_clusters <- params$characterized_clusters_strategy(robust_clusters, data.iteration, params)

  if (length(characterized_clusters) == 0) {return(list())}
  clusters <- add_leftover_cluster(population, characterized_clusters, data.iteration, params)

  if (length(clusters) == 2) {
    characterized_clusters <- params$characterized_clusters_strategy(clusters, data.iteration, params)
    if (length(characterized_clusters) < 2) {return(list())}}

  if (figures) {
    plot <- draw_marker_genes(clusters)
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_marker_genes.pdf"))
    print(plot)
    grDevices::dev.off()}
  return(clusters)
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

draw_marker_genes <- function(characterized_clusters) {
  #' Get an upset-plot representing the marker genes predicted in each meta-cluster.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' @return a plot.
  #'
  #' @import ggVennDiagram
  #' @import SCpubr
  #'
  #' @export
  #'
  get_markers.cluster <- function(cluster) {names(cluster$markers)}
  markers <- sapply(X=characterized_clusters, FUN=get_markers.cluster)
  labels <- sapply(X=characterized_clusters, FUN="[[", "label")
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
