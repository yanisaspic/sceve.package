"Main functions called to condut a scEVE clustering analysis.

	2025/01/16 @yanisaspic"

get_default_parameters <- function() {
  #' Get the default parameters of the scEVE algorithm. These parameters are:
  #' `random_state` the integer seed used to have deterministic results.
  #' `figures_path` & `sheets_path` the default paths where figures and result sheets are stored, respectively.
  #'
  #' `selected_genes_strategy` a function called to select a limited pool of genes for a clustering iteration.
  #' This function is called by `extract_data()`, and it expects two positional arguments: `expression` and `params`, respectively.
  #' By default, the n most variable genes for a clustering iteration are selected, with n=500.
  #' It outputs a vector of genes.
  #'
  #' `base_clusters_strategy` a function called to predict base clusters with multiple clustering methods.
  #' This function is called by `get_base_clusters()`, and it expects two positional arguments: `data.iteration` and `params`, respectively.
  #' By default, four clustering methods are used: densityCut, monocle3, Seurat and SHARP.
  #' It outputs a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' `robust_clusters_strategy` a function called to identify robust clusters from the similarity subgraphs.
  #' This function is called by `get_meta_clusters()`, and it expects two positional arguments: `subgraphs` and `params`, respectively.
  #' By default, a robustness threshold is used to filter out subgraphs, with threshold=0.33.
  #' It outputs a list where every element is a pool of cells grouped together by multiple clustering methods.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label` and `robustness`.
  #'
  #' `marker_genes_strategy` a function called to predict the marker genes of every meta-cluster.
  #' This function is called by `get_characterized_clusters()`, and it expects three positional arguments: `meta_clusters`, `data.iteration` and `params`.
  #' By default, an over-representation test is conducted on every cluster to identify genes that are
  #' frequently expressed in it, but rarely expressed in the complete pool of cells.
  #' It outputs a list where every element is a pool of cells.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers` and `robustness`.
  #'
  #' `characterized_clusters_strategy` a function called to identify characterized clusters.
  #' This function is called by `get_characterized_clusters()`, and it expects two positional arguments: `meta_clusters` and `params`.
  #' By default, a characterized cluster is a meta-cluster with at least 23 specific marker genes.
  #' It outputs a list where every element is a pool of cells.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `cells`, `clustering_methods`, `label`, `markers`, `robustness` and `specific_markers`.
  #'
  #' `leftover_cells_strategy` a function called to manage the cells in the leftover cluster of a clustering iteration.
  #' This function is called by ``
  #'
  #' @return a list of parameters.
  #'
  #' @export
  #'
  params <- list(
    random_state=0,
    figures_path="./scEVE/figures",
    sheets_path="./scEVE/sheets",
    selected_genes_strategy=get_selected_genes.n_most_variable,
    base_clusters_strategy=get_base_clusters.default_methods,
    robust_clusters_strategy=get_robust_clusters.robustness_threshold,
    marker_genes_strategy=get_marker_genes.over_representation,
    characterized_clusters_strategy=get_characterized_clusters.markers_threshold,
    leftovers_strategy="default",
  )
  return(params)
}

initialize_records <- function(expression.init) {
  #' Get a named list, with three data.frames: `cells`, `markers` and `meta`.
  #' - `cells` associates cells to their predicted populations.
  #' Its rows are cells, its columns are predicted populations, and cell memberships are reported in the table.
  #' - `markers` associates predicted populations to their marker genes.
  #' Its rows are genes, its columns are predicted populations, and gene representations are reported in the table.
  #' - `meta` associates predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a named list, with three data.frames: `cells`, `markers` and `meta`.
  #'
  #' @export
  #'
  cells <- data.frame(
    C=as.numeric(rep(1, ncol(expression.init))),
    row.names=colnames(expression.init)
  )
  markers <- data.frame(
    C=as.numeric(rep(0, nrow(expression.init))),
    row.names=rownames(expression.init)
  )
  meta <- data.frame(
    size=as.numeric(ncol(expression.init)),
    robustness=0,
    parent=NA,
    clustering_status="PENDING",
    row.names="C"
  )
  records <- list(cells=cells, markers=markers, meta=meta)
  return(records)
}

get_SeuratObject.init <- function(expression.init, params) {
  #' Get a SeuratObject from a scRNA-seq dataset of raw count expression, without selected genes.
  #' This function is used once prior to a scEVE clustering analysis in order to draw the
  #' extracted data at each clustering iteration (cf. `sceve::draw_extracted_data()`).
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `sceve::get_default_parameters()`).
  #'
  #' @return a SeuratObject, on which the function RunUMAP() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  #' @export
  #'
  SeuratObject.init <- Seurat::CreateSeuratObject(expression.init)
  selected_genes <- params$selected_genes_strategy(expression.init)
  Seurat::VariableFeatures(SeuratObject.init) <- selected_genes
  SeuratObject.init <- Seurat::NormalizeData(SeuratObject.init)
  SeuratObject.init <- Seurat::ScaleData(SeuratObject.init,
                                         features=Seurat::VariableFeatures(SeuratObject.init))
  SeuratObject.init <- Seurat::RunUMAP(SeuratObject.init,
                                       features=Seurat::VariableFeatures(SeuratObject.init),
                                       seed.use=params$random_state)
  return(SeuratObject.init)
}

#' #'
#' #' scEVE.iteration <- function(population, expression.init, SeurObj.init, records, params,
#' #'                             figures, sheets) {
#' #'   #' Attempt to cluster a specific cell population using the scEVE algorithm.
#' #'   #'
#' #'   #' @param population: a character. It corresponds to the cell population that scEVE will attempt to cluster.
#' #'   #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes.
#' #'   #' Its rows are genes and its columns are cells.
#' #'   #' @param SeurObj.init: a SeuratObject generated from expression.init.
#' #'   #' A U-MAP is already applied on this object.
#' #'   #' @param records: a named list, with three data.frames: `cells`, `markers` and `meta`.
#' #'   #' @param params: a list of parameters. See sceve::get_default_hyperparameters().
#' #'   #' @param figures: a boolean that indicates if figures should be drawn to explain the clustering iteration.
#' #'   #' @param sheets: a boolean that indicates if the results of the clustering iteration should be saved in Excel sheets.
#' #'   #'
#' #'   #' @return a named list, with three data.frames: `cells`, `markers` and `meta`.
#' #'   #'
#' #'
#' #'   while (TRUE) {
#' #'     data.loop <- trim_data(population, expression.init, SeurObj.init, records, params,
#' #'                            figures)
#' #'     if (length(data.loop)==0){break()}
#' #'     # _______________________________________if too little cells, do not try to cluster
#' #'
#' #'     clusterings <- get_clusterings(population, expression.init, data.loop,
#' #'                                    params, figures)
#' #'     seeds <- get_seeds(population, expression.init, data.loop, clusterings, records,
#' #'                        params, figures)
#' #'     if (length(seeds)==0){break()}
#' #'     # ______________________________if too little consensus for every seed, do not characterize
#' #'
#' #'     data.loop$occurrences.loop <- get_occurrences(data.loop$ranked_genes.loop)
#' #'     seeds <- get_genes(population, data.loop, seeds, params, figures)
#' #'     if (length(seeds)==0){break()}
#' #'     # _____________________________if too little characterization for any seed, do not report
#' #'
#' #'     records <- update_records(population, data.loop, seeds, records, params)
#' #'     if (sheets) {write.xlsx(records, params$sheets_path, rowNames=TRUE)}
#' #'     break()
#' #'   }
#' #'
#' #'   if (figures){merge_pdfs(population, params)}
#' #'   return(records)
#' #' }
#' #'
#' #' do_scEVE <- function(expression.init, params=get_default_hyperparameters(),
#' #'                      figures=TRUE, sheets=TRUE) {
#' #'   #' Attempt to cluster a scRNA-seq dataset using the scEVE algorithm.
#' #'   #'
#' #'   #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes.
#' #'   #' Its rows are genes and its cols are cells.
#' #'   #' @param params: a list of parameters (cf. `sceve::get_default_parameters()`).
#' #'   #' @param figures: a boolean that indicates if figures should be drawn to explain the clustering iterations.
#' #'   #' @param sheets: a boolean that indicates if the results of the clustering iteration should be saved in Excel sheets.
#' #'   #'
#' #'   #' @return a named list, with two elements:
#' #'   #' - `records`: a named list, with three data.frames: `cells`, `markers` and `meta`.
#' #'   #' - `labels`: a named factor. Its names are cells and its values are population labels.
#' #'   #'
#' #'
#' #'   #_________________________________________________________________________init
#' #'   records <- initialize_records(expression.init)
#' #'   if (figures) {
#' #'     dir.create(params$figures_path)
#' #'     SeuratObject.init <- get_SeuratObject.init(expression.init, params)
#' #'   } else {SeuratObject.init <- NA}
#' #'   population <- "C"
#' #'
#' #'   #____________________________________________________________________main loop
#' #'   while (!is.na(population)) {
#' #'     records <- scEVE.iteration(population, expression.init, SeurObj.init, records, params,
#' #'                                figures, sheets)
#' #'     records$meta[population, "to_dig"] <- FALSE
#' #'     population <- get_undug_population(records)
#' #'     gc()
#' #'   }
#' #'
#' #'   #_______________________________________________________________________finale
#' #'   is_marker <- function(row) {sum(row) > 0}
#' #'   records$markers <- records$markers[apply(X=records$markers, MARGIN=1, FUN=is_marker),]
#' #'   if (sheets) {write.xlsx(records, params$sheets_path, rowNames=TRUE)}
#' #'   results <- list(records=records, preds=factor(get_leaves(records$cells)))
#' #'   return(results)
#' #' }
