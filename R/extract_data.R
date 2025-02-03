"Functions called to set-up the data of a scEVE clustering iteration.

	2025/01/16 @yanisaspic"

get_selected_genes.n_most_variable <- function(expression, params, n_genes=500) {
  #' Get the n most variable genes in a scRNA-seq dataset.
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
  expression.selected <- expression.selected[, colSums(expression.selected)>0]
    # ignore cells that do not express any selected genes.
  return(expression.selected)
}

ignore_dropped_genes <- function(ranking_of_genes, expression.selected) {
  #' Remove genes that are not expressed in a cell, from a table ranking gene expression.
  #'
  #' @param ranking_of_genes a data.frame associating cells to their genes, ranked by descending order of expression.
  #' Its rows are ranks, its columns are cells, and genes are reported in the table.
  #' @param expression.selected a scRNA-seq dataset of raw count expression, with selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a data.frame associating cells to their genes, ranked by descending order of expression.
  #' Its rows are ranks, its columns are cells, and either genes or NA values are reported in the table.
  #'
  n_genes <- nrow(ranking_of_genes)
  get_n_expressed_genes.cell <- function(expr.cell){sum(expr.cell>0)}
  n_expressed_genes <- apply(X=expression.selected, MARGIN=2, FUN=get_n_expressed_genes.cell)
  for (cell in colnames(ranking_of_genes)) {
    dropped_genes <- (n_expressed_genes[cell]+1):n_genes
    ranking_of_genes[dropped_genes, cell] <- NA
  }
  return(ranking_of_genes)
}

get_SeuratObject.selected <- function(expression.selected) {
  #' Get a SeuratObject from a scRNA-seq dataset of raw count expression, with selected genes.
  #' This function is used as a pre-processing step for the Seurat and monocle3 clustering algorithms.
  #'
  #' @param expression.selected a scRNA-seq dataset of raw count expression, with selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a SeuratObject, on which the function ScaleData() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  SeuratObject <- Seurat::CreateSeuratObject(expression.selected)
  Seurat::VariableFeatures(SeuratObject) <- rownames(expression.selected)
  SeuratObject <- Seurat::NormalizeData(SeuratObject)
  SeuratObject <- Seurat::ScaleData(SeuratObject,
                                    features=Seurat::VariableFeatures(SeuratObject))
  return(SeuratObject)
}

get_ranking_of_genes <- function(expression.selected) {
  #' Rank the genes of every cell by descending order of expression, and get the resulting table.
  #' Only expressed genes are reported, and the number of genes expressed varies across cells.
  #' Hence, some table cells include NA values.
  #'
  #' @param expression.selected a scRNA-seq dataset of raw count expression, with selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a data.frame associating cells to their genes, ranked by descending order of expression.
  #' Its rows are ranks, its columns are cells, and genes are reported in the table.
  #'
  get_ranks.cell <- function(expr.cell){rank(-expr.cell, ties.method="random")}
  get_genes.cell <- function(ranks.cell){names(sort(ranks.cell))}
  get_ranking_of_genes.cell <- function(expr.cell){get_genes.cell(get_ranks.cell(expr.cell))}
  ranking_of_genes <- apply(X=expression.selected, MARGIN=2, FUN=get_ranking_of_genes.cell)
  ranking_of_genes <- ignore_dropped_genes(ranking_of_genes, expression.selected)
  return(ranking_of_genes)
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
  #' @return a named list, with four names: `expression`, `SeuratObject`, `expression.init` and `ranking_of_genes.init`.
  #' The two first elements correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject. and the ranking of genes generated from this matrix.
  #' The two last elements correspond to the full scRNA-seq expression matrix and the ranking of genes generated from this matrix.
  #'
  #' @import Seurat
  #' @import glue
  #' @import grDevices
  #' @import qpdf
  #'
  #' @export
  #'
  cells_of_population <- get_cells_of_population(population, records$cells)
  if (length(cells_of_population)>=100) {
    # the parameter sncells of the clustering algorithm SHARP is set to 100 by default.
    expression.population <- expression.init[, cells_of_population]
    expression.population <- get_expression.selected(expression.population, params)
    SeuratObject.population <- get_SeuratObject.selected(expression.population)
    SeuratObject.population <- Seurat::RunUMAP(SeuratObject.population,
                                               features=Seurat::VariableFeatures(SeuratObject.population),
                                               seed.use=params$random_state)
    ranking_of_genes.init <- get_ranking_of_genes(expression.init[rownames(expression.population), ])
    data.iteration <- list(expression=expression.population,
                           SeuratObject=SeuratObject.population,
                           expression.init=expression.init,
                           ranking_of_genes.init=ranking_of_genes.init)}
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
