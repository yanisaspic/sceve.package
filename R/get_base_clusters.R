"Functions used to get base clusters.

	2025/01/21 @yanisaspic"

use_Seurat <- function(SeuratObject.population, params) {
  #' Predict clusters with the Seurat method.
  #'
  #' @param SeuratObject.population a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import Seurat
  #'
  if (!"umap" %in% names(SeuratObject.population@reductions)) {
  SeuratObject.population <- Seurat::RunPCA(SeuratObject.population,
                                            features=Seurat::VariableFeatures(SeuratObject.population),
                                            seed.use=params$random_state)}
  SeuratObject.population <- Seurat::FindNeighbors(SeuratObject.population,
                                                   features=Seurat::VariableFeatures(SeuratObject.population))
  SeuratObject.population <- Seurat::FindClusters(SeuratObject.population,
                                                  random.seed=params$random_state)
  predictions <- Seurat::Idents(SeuratObject.population)
  return(predictions)
}

use_monocle3 <- function(SeuratObject.population, params) {
  #' Predict clusters with the monocle3 method.
  #'
  #' @param SeuratObject.population a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import monocle3
  #' @import R.utils
  #' @import Seurat
  #' @import SeuratWrappers
  #'
  if (!"umap" %in% names(SeuratObject.population@reductions)) {
    SeuratObject.population <- Seurat::RunUMAP(SeuratObject.population,
                                               features=Seurat::VariableFeatures(SeuratObject.population),
                                               seed.use=params$random_state)}
  CDSObject <- SeuratWrappers::as.cell_data_set(SeuratObject.population)
  CDSObject <- monocle3::cluster_cells(CDSObject, random_seed=params$random_state)
  predictions <- CDSObject@clusters@listData$UMAP$clusters
  return(predictions)
}

get_formatted_predictions <- function(cells, predictions) {
  #' Format the cluster predictions of a method to obtain a standard output.
  #'
  #' @param cells a vector of cells names.
  #' @param predictions a vector of cluster predictions.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import stats
  #'
  predictions <- stats::setNames(predictions, cells)
  predictions <- factor(predictions)
  return(predictions)
}

use_SHARP <- function(expression.population, params) {
  #' Predict clusters with the SHARP method.
  #'
  #' @param expression.population a scRNA-seq dataset of raw count expression, with selected genes.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import clues
  #' @import SHARP
  #'
  results <- SHARP::SHARP(scExp=expression.population, exp.type="count",
                          n.cores = 2, rN.seed=params$random_state)
  predictions <- get_formatted_predictions(cells=colnames(expression.population),
                                           predictions=results$pred_clusters)
  return(predictions)
}

use_densityCut <- function(logtpm.population, params) {
  #' Predict clusters with the densityCut method.
  #'
  #' @param logtpm.population a scRNA-seq dataset of log-transformed transcripts per millions (tpm),
  #' with selected genes.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import densitycut
  #'
  set.seed(params$random_state)
  data <- t(logtpm.population) # densityCut expects cells as rows and genes as columns.
  results <- densitycut::DensityCut(t(logtpm.population), show.plot = FALSE)
  predictions <- get_formatted_predictions(cells=colnames(logtpm.population), predictions=results$cluster)
  return(predictions)
}

get_data_input <- function(method, data) {
  #' Get the data input relevant to a clustering method.
  #'
  #' @param method a clustering method used by scEVE.
  #' Valid methods include: densityCut, monocle3, Seurat and SHARP.
  #' @param data a named list, with three names: `expression`, `SeuratObject` and `logtpm`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population, its SeuratObject and its log2-transformed TPM values.
  #'
  #' @return one of `data$expression`, `data$SeuratObject` or `data$logtpm`.
  #'
  data_inputs <- list(Seurat=data$SeuratObject,
                      monocle3=data$SeuratObject,
                      SHARP=data$expression,
                      densityCut=data$logtpm)
  input <- data_inputs[[method]]
  return(input)
}

get_base_clusters.default_methods <- function(data.iteration, params,
                                              clustering_methods=c("densityCut", "monocle3", "Seurat", "SHARP")) {
  #' Apply multiple clustering methods on scRNA-seq data to predict base clusters.
  #'
  #' @param data.iteration a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param clustering_methods a vector of valid clustering methods.
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @import glue
  #' @import scater
  #'
  #' @export
  #'
  logtpm.population <- log2(scater::calculateTPM(data.iteration$expression) + 1)
  data.iteration[["logtpm"]] <- logtpm.population

  get_predictions.method <- function(method) {
    f <- get(glue("use_{method}"))
    predictions <- f(get_data_input(method, data.iteration), params)
    return(predictions)
  }
  base_clusters <- sapply(X=clustering_methods, FUN=get_predictions.method)
  gc()

  rename_clusters.method <- function(method) {
    rename.prediction <- function(prediction) {glue("{method}_{prediction}")}
    column <- sapply(X=base_clusters[, method], FUN=rename.prediction)
    return(column)
  }
  base_clusters <- sapply(X=clustering_methods, FUN=rename_clusters.method)
  return(base_clusters)
}

draw_base_clusters <- function(SeuratObject.population, base_clusters) {
  #' Get composite U-MAP plots representing the base clusters predicted by each clustering method.
  #'
  #' @param SeuratObject.population a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param base_clusters a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
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
  get_number <- function(renamed_cluster){
    number <- strsplit(renamed_cluster, split="_")[[1]][2]
  }
  base_clusters <- apply(base_clusters, c(1,2), get_number)

  get_plot.method <- function(method) {
    cells <- rownames(base_clusters)
    numbers <- factor(base_clusters[, method])
    SeuratObject.population@active.ident <- stats::setNames(numbers, cells)
    n_clusters <- length(levels(numbers))
    plot <- SCpubr::do_DimPlot(SeuratObject.population) +
      ggplot2::ggtitle(glue("{method}: {n_clusters} cluster(s)")) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title=element_text(hjust=0.5, margin=margin(1, 0, 0, 0)),
                     panel.background=element_rect(fill="lightgrey"), legend.position="none",
                     axis.title=element_blank())
    return(plot)
  }

  clustering_methods <- colnames(base_clusters)
  plots <- lapply(X=clustering_methods, FUN=get_plot.method)
  names(plots) <- clustering_methods
  composite_plot <- do.call(gridExtra::grid.arrange, plots)
  return(composite_plot)
}

get_base_clusters <- function(population, data.iteration, params, figures) {
  #' Get base clusters predicted with multiple clustering methods.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param data.iteration a named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @import glue
  #' @import grDevices
  #'
  #' @export
  #'
  base_clusters <- params$base_clusters_strategy(data.iteration, params)
  if (figures) {
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_base_clusters.pdf"))
    plot <- draw_base_clusters(data.iteration$SeuratObject, base_clusters)
    print(plot)
    grDevices::dev.off()
  }
  return(base_clusters)
}
