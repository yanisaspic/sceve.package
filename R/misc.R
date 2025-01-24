"Miscellaneous functions called multiple times in the scEVE algorithm.

	2025/01/16 @yanisaspic"

get_resolution <- function(population) {
  #' Get the resolution level of a predicted population.
  #'
  #' @param population a character. It corresponds to a cell population predicted by scEVE (e.g. C.2).
  #'
  #' @return an integer.
  #'
  resolution <- length(strsplit(population, split=".", fixed=TRUE)[[1]])
  return(resolution)
}

get_populations_at_resolution <- function(resolution, records.cells) {
  #' Get all the populations at a specific resolution.
  #' The root population is resolution 1, and its children populations are resolution 2, etc.
  #'
  #' @param resolution an integer.
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of population labels.
  #'
  populations <- colnames(records.cells)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  populations_at_resolution <- populations[resolutions==resolution]
  return(populations_at_resolution)
}

get_max_resolution <- function(records.cells) {
  #' Get the maximum clustering resolution attained in a scEVE clustering analysis.
  #'
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return an integer.
  #'
  populations <- colnames(records.cells)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  max_resolution <- max(resolutions)
  return(max_resolution)
}

get_cells_of_population <- function(population, records.cells) {
  #' Get every cell belonging in a predicted population.
  #'
  #' @param population a character. It corresponds to a cell population predicted by scEVE (e.g. C.2).
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of cells.
  #'

  if (population=="C") {return(rownames(records.cells))}
  # all cells belong to the root population

  cell_is_in_population <- function(cell_membership) {
    max_membership <- max(cell_membership)
    main_population <- names(which.max(cell_membership))
    is_in_population <- (main_population == population) &
      (max_membership > 0)
  }

  resolution <- get_resolution(population)
  records.cells <- records.cells[, get_populations_at_resolution(resolution, records.cells)]
  cells_are_in_population <- apply(X=records.cells, MARGIN=1, FUN=cell_is_in_population)
  cells_of_population <- rownames(records.cells[cells_are_in_population, ])
  return(cells_of_population)
}

get_gene_occurrences <- function(ranking_of_genes) {
  #' Get the gene occurrences with regards to the gene sampling effort.
  #' The occurrences correspond to the number of cells expressing a gene, and the
  #' sampling effort corresponds to the number of genes sampled in a cell.
  #' The genes are sampled in a descending order, so that the n most expressed genes
  #' are sampled with an effort of n.
  #' e.g.: only the most expressed gene is sampled with a sampling effort of 1,
  #' whereas the 3 most expressed genes are sampled with a sampling effort of 3.
  #'
  #' @param ranking_of_genes a data.frame associating cells to their genes, ranked by descending order of expression.
  #' Its rows are ranks, its columns are cells, and genes are reported in the table.
  #'
  #' @return a data.frame associating gene occurrences to the sampling effort.
  #' Its rows are genes, its columns are sampling efforts, and gene occurrences are reported in the table.
  #'
  max_effort <- max(apply(X=!is.na(ranking_of_genes), MARGIN=2, FUN=sum))
  genes <- sort(unique(as.vector(ranking_of_genes)))

  get_occurrences_at_effort <- function(effort) {
    data <- ranking_of_genes[1:effort,]
    occurrences_at_effort <- table(as.vector(data))
    missing_genes <- setdiff(genes, names(occurrences_at_effort))
    occurrences_at_effort[missing_genes] <- 0
    return(occurrences_at_effort[genes])
  }

  results <- lapply(X=1:max_effort, FUN=get_occurrences_at_effort)
  occurrences <- as.data.frame(do.call(cbind, results))
  occurrences <- as.data.frame(occurrences)
  colnames(occurrences) <- 1:max_effort
  occurrences[] <- lapply(occurrences, as.numeric)

  return(occurrences)
}

get_records <- function(sheets_path) {
  #' Load the records of a scEVE clustering analysis.
  #'
  #' @param sheets_path a path where result sheets are stored.
  #'
  #' @return a named list, with three data.frames: `cells`, `markers` and `meta`.
  #'
  #' @import openxlsx
  #'
  sheets_names <- openxlsx::getSheetNames(sheets_path)
  get_sheet <- function(sheets_name) {openxlsx::read.xlsx(sheets_path, sheet=sheets_name, rowNames=TRUE)}
  records <- sapply(X=sheets_names, FUN=get_sheet)
  return(records)
}
