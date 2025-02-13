"Functions called to set-up the data of a scEVE clustering iteration.

	2025/02/ @yanisaspic"

get_selected_genes.n_most_variable <- function(expression, params, n_genes=1000) {
  #' Get the n most variable genes in a scRNA-seq dataset.
  #'
  #' The variable genes are identified by calling the function FindVariableFeatures() of Seurat.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param n_genes the number of highly variable genes to sample.
  #'
  #' @return a vector of genes.
  #'
  #' @import Seurat
  #'
  #' @export
  #'
  SeurObj <- Seurat::CreateSeuratObject(expression)
  SeurObj <- Seurat::FindVariableFeatures(SeurObj, nfeatures = n_genes)
  return(Seurat::VariableFeatures(SeurObj))
}

get_expression.selected <- function(expression, params) {
  #' Get a scRNA-seq dataset of raw count expression, with selected genes.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a scRNA-seq dataset of raw count expression, with selected genes.
  #'
  selected_genes <- params$selected_genes_strategy(expression, params)
  expression.selected <- expression[selected_genes,]
  return(expression.selected)
}

get_SeuratObject.selected <- function(expression.selected, params) {
  #' Get a SeuratObject from a scRNA-seq dataset of raw count expression, with selected genes.
  #'
  #' This function is used as a pre-processing step for the Seurat and monocle3 clustering algorithms.
  #'
  #' @param expression.selected a scRNA-seq dataset of raw count expression, with selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a SeuratObject, on which the function ScaleData() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  SeurObj <- Seurat::CreateSeuratObject(expression.selected)
  Seurat::VariableFeatures(SeurObj) <- rownames(expression.selected)
  SeurObj <- Seurat::NormalizeData(SeurObj)
  SeurObj <- Seurat::ScaleData(SeurObj, features=Seurat::VariableFeatures(SeurObj))
  SeurObj <- Seurat::RunPCA(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=params$random_state)
  SeurObj <- Seurat::RunUMAP(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=params$random_state)
  return(SeurObj)
}

extract_data <- function(population, expression.init, SeuratObject.init, records, params, figures) {
  #' Extract a data subset with a specific cell population and their most variableg genes.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param SeuratObject.init a SeuratObject generated from expression.init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #'
  #' @return a named list, with two names: `expression` and `SeuratObject`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.
  #'
  #' @import Seurat
  #' @import glue
  #' @import grDevices
  #' @import qpdf
  #'
  #' @export
  #'
  cells_of_population <- get_cells_of_population(population, records$cells)
  if (length(cells_of_population) >= params$minimum_cells) {
    # the parameter sncells of the clustering algorithm SHARP is set to 100 by default.
    expression.population <- expression.init[, cells_of_population]
    expression.population <- get_expression.selected(expression.population, params)
    SeuratObject.population <- get_SeuratObject.selected(expression.population, params)
    data.iteration <- list(expression=expression.population,
                           SeuratObject=SeuratObject.population)}
  else {data.iteration <- list()}

  if (figures) {
    plot <- draw_extracted_data(population, SeuratObject.init, records)
    grDevices::pdf(file = glue::glue("{params$figures_path}/{population}_extract_data.pdf"))
    print(plot)
    grDevices::dev.off()
  }
  return(data.iteration)
}

draw_extracted_data <- function(population, SeuratObject.init, records) {
  #' Get a U-MAP plot representing the pool of cells used in the clustering iteration.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param SeuratObject.init a SeuratObject generated from expression.init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a plot.
  #'
  #' @import viridis
  #' @import forcats
  #' @import assertthat
  #' @import ggplotify
  #' @import SCpubr
  #' @import ggplot2
  #'
  #' @export
  #'
  cells_of_population <- get_cells_of_population(population, records$cells)
  plot <- SCpubr::do_DimPlot(SeuratObject.init, cells.highlight=cells_of_population) +
    ggplot2::ggtitle(population) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.background=element_rect(fill="lightgrey"),
                   axis.title=element_blank(), legend.position="bottom")
  return(plot)
}
