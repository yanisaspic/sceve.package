"Main functions called to condut a scEVE clustering analysis.

	2025/01/16 @yanisaspic"

get_default_parameters <- function() {
  #' Get the default parameters of the scEVE algorithm. These parameters are:
  #' `random_state` the integer seed used to have deterministic results.
  #' `robustness_threshold` the threshold used to identify robust clusters.
  #' `minimum_cells` the minimum number of cells expected in a population before attempting to sub-cluster it.
  #' `figures_path` & `sheets_path` the default paths where figures and result sheets are stored, respectively.
  #' `selected_genes_strategy` a function called to select a limited pool of genes for a clustering iteration.
  #' `base_clusters_strategy` a function called to predict base clusters with multiple clustering methods.
  #' `marker_genes_strategy` a function called to predict the marker genes of every meta-cluster.
  #' `characterized_clusters_strategy` a function called to identify characterized clusters.
  #' `cluster_memberships_strategy` a function called to quantify the
  #' Detailed information regarding the strategy parameters are available in the vignette of the package.
  #'
  #' @return a list of parameters.
  #'
  #' @export
  #'
  params <- list(random_state=1, robustness_threshold=0.33, minimum_cells=100,
                 figures_path="./scEVE", sheets_path="./scEVE/records.xlsx",
                 selected_genes_strategy=get_selected_genes.n_most_variable,
                 base_clusters_strategy=get_base_clusters.default_methods,
                 marker_genes_strategy=add_marker_genes.over_representation,
                 characterized_clusters_strategy=get_characterized_clusters.markers_threshold,
                 cluster_memberships_strategy=get_cluster_memberships.binary_membership)
  return(params)
}

initialize_records <- function(expression.init) {
  #' Get a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' - `cells` associates cells to their predicted populations.
  #' Its rows are cells, its columns are predicted populations, and cell memberships are reported in the table.
  #' - `markers` associates predicted populations to their marker genes.
  #' Its rows are genes, its columns are predicted populations, and characterization powers are reported in the table.
  #' - `meta` associates predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #' - `methods` associates predicted populations to the clustering methods leveraged to predict them.
  #' Its rows are clustering methods, its columns are predicted populations, and binary values are reported in the table.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @export
  #'
  cells <- data.frame(C=as.numeric(rep(1, ncol(expression.init))), row.names=colnames(expression.init))
  markers <- data.frame(C=as.numeric(rep(0, nrow(expression.init))), row.names=rownames(expression.init))
  meta <- data.frame(size=as.numeric(ncol(expression.init)), robustness=0, parent=NA,
                     clustering_status="PENDING", row.names="C")
  methods <- data.frame()
  records <- list(cells=cells, markers=markers, meta=meta, methods=methods)
  return(records)
}

get_SeuratObject.init <- function(expression.init) {
  #' Get a SeuratObject from a scRNA-seq dataset of raw count expression, without selected genes.
  #' This function is used once prior to a scEVE clustering analysis in order to draw the
  #' extracted data at each clustering iteration (cf. `sceve::draw_extracted_data()`).
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a SeuratObject, on which the function RunUMAP() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  #' @export
  #'
  SeuratObject.init <- Seurat::CreateSeuratObject(expression.init)
  SeuratObject.init <- Seurat::FindVariableFeatures(SeuratObject.init, nfeatures = 1000)
  SeuratObject.init <- Seurat::NormalizeData(SeuratObject.init)
  SeuratObject.init <- Seurat::ScaleData(SeuratObject.init,
                                         features=Seurat::VariableFeatures(SeuratObject.init))
  SeuratObject.init <- Seurat::RunUMAP(SeuratObject.init,
                                       features=Seurat::VariableFeatures(SeuratObject.init),
                                       seed.use=1)
  return(SeuratObject.init)
}

get_pending_population <- function(records) {
  #' Get a cell population for which no scEVE clustering iteration has been attempted.
  #'
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a character.
  #'
  meta <- records$meta
  pending_populations <- rownames(meta[meta$clustering_status=="PENDING", ])
  population <- pending_populations[1]
  return(population)
}

sceve.iteration <- function(population, expression.init, SeuratObject.init, records, params, figures, sheets) {
  #' Attempt to cluster a specific cell population using the scEVE algorithm.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param SeuratObject.init a SeuratObject generated from expression.init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #' @param sheets a boolean that indicates if the results of the clustering iteration should be saved in Excel sheets.
  #'
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @import openxlsx
  #'
  while(TRUE) {
    data.iteration <- extract_data(population, expression.init, SeuratObject.init, records, params, figures)
    if (length(data.iteration) == 0) {break()} # the cell population is too small
    base_clusters <- get_base_clusters(population, data.iteration, params, figures)
    meta_clusters <- get_meta_clusters(population, base_clusters, data.iteration, records, params, figures)
    if (length(meta_clusters) == 1) {break()} # multiple sub-clusters are not predicted
    characterized_clusters <- get_characterized_clusters(population, meta_clusters, data.iteration, params, figures)
    if (length(characterized_clusters) == 0) {break()}  # the sub-clusters are homogenous
    records <- report_iteration(population, characterized_clusters, data.iteration, records, params)
    break()}
  records$meta[population, "clustering_status"] <- "COMPLETE"
  if (sheets) {openxlsx::write.xlsx(records, params$sheets_path, rowNames=TRUE)}
  if (figures) {merge_drawings(population, params)}
  return(records)
}

sceve <- function(expression.init, params=get_default_parameters(), figures=TRUE, sheets=TRUE) {
  #' Conduct a clustering analysis with the scEVE framework.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering iteration.
  #' @param sheets a boolean that indicates if the results of the clustering iteration should be saved in Excel sheets.
  #'
  #' @return a named list, with two elements: `records` and `preds`.
  #' `records` is a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #' `preds` is a named factor associating cells to their predicted clusters.
  #'
  #' @import openxlsx
  #'
  #' @export
  #'
  records <- initialize_records(expression.init)
  if (figures) {
    dir.create(params$figures_path)
    SeuratObject.init <- get_SeuratObject.init(expression.init)}
  else {SeuratObject.init <- NA}
  population <- "C"

  while (!is.na(population)) {
    records <- sceve.iteration(population, expression.init, SeuratObject.init,
                               records, params, figures, sheets)
    population <- get_pending_population(records)}

  gene_is_marker <- function(gene) {sum(gene) > 0}
  records$markers <- records$markers[apply(X=records$markers, MARGIN=1, FUN=gene_is_marker),]
  if (sheets) {openxlsx::write.xlsx(records, params$sheets_path, rowNames=TRUE)}
  results <- list(records=records, preds=factor(get_leaf_clusters(records$cells)))
  return(results)
}
