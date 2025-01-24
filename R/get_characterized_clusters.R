"Functions called to identify characterized clusters, i.e. clusters with sufficient marker genes.

	2025/01/24 @yanisaspic"

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggplotify)
  library(ggVennDiagram)
  library(glue)
  library(stats)
})
source("./src/scEVE/utils/misc.R")

add_occurrences_to_seeds <- function(ranked_genes, seeds) {
  #' Add the occurrences respective to each seed. This information is used for the overrepresentation test.
  #'
  #' @param ranked_genes: a data.frame where: ranks are rows | cells are cols | cells are genes.
  #' @param seeds: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #'
  #' @return a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'genes'.
  #'
  get_occurrences.seed <- function(seed){get_occurrences(ranked_genes[, seed$cells])}
  occurrences <- lapply(X=seeds, FUN=get_occurrences.seed)
  for (i in 1:length(seeds)) {
    seeds[[i]]$occurrences <- occurrences[[i]]
  }
  return(seeds)
}

get_efforts.plot <- function(efforts.frame) {
  #' Get a boxplot corresponding to the effort w.r.t. the group of cells.
  #' The whiskers correspond to the minimum and maximum values.
  #'
  #' @param efforts.frame: a data.frame with two columns: 'effort' and 'seed'.
  #'
  #' @return a boxplot where x=seed, y=effort and group=seed.
  #'
  efforts.plot <- ggplot(data=efforts.frame,
                         aes(x=seed, y=effort, group=seed, fill=seed)) +
    geom_boxplot(coef=NULL) +
    ggtitle("Sampling effort of cells, per consensus cluster.")
  return(efforts.plot)
}

test_overrepresentation <- function(q, m, N, k) {
  #' Get the probability of observing x successes or more in a sample.
  #'
  #' @param q: number of successes in sample.
  #' @param m: number of successes in reference.
  #' @param n: size of the reference.
  #' @param k: size of the sample.
  #'
  #' @return a numeric.
  #'

  # phyper tests Pr(X > q): so use q-1 #
  ######################################
  n <- N - m
  pvalue <- phyper(q-1, m, n, k, lower.tail=FALSE)
  return(pvalue)
}

test_overrepresentation.gene <- function(gene, overrepresentation.frame) {
  #' Test if a gene is over-represented in a sample w.r.t. a reference.
  #'
  #' @param gene: a character.
  #' @param overrepresentation.frame: a data.frame with two columns: 'sample' and 'population'. The row names are genes.
  #'
  #' @return a numeric.
  #'
  q <- overrepresentation.frame[gene, "sample"]
  m <- overrepresentation.frame[gene, "population"]
  N <- sum(overrepresentation.frame$population)
  k <- sum(overrepresentation.frame$sample)
  pvalue <- test_overrepresentation(q, m, N, k)
  return(pvalue)
}

test_overrepresentation.seed <- function(seed, occurrences.population) {
  #' Test which genes in a seed are overrepresented w.r.t. minimal effort in the seed.
  #'
  #' @param seed: a list with four keys: 'consensus', 'cells', 'clusters' and 'genes'.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #'
  #' @return a data.frame with four columns: 'sample', 'population', 'pvalue' and 'adj_pvalue'
  #'
  occurrences.seed <- seed$occurrences
  occurrences.population <- occurrences.population[rownames(occurrences.seed),]
  overrepresentation.frame <- data.frame(sample=occurrences.seed[, length(occurrences.seed)],
                                         population=occurrences.population[, length(occurrences.population)])
  rownames(overrepresentation.frame) <- rownames(occurrences.seed)

  ### p-values w/ correction ###
  ##############################
  pvalues <- sapply(X=rownames(overrepresentation.frame),
                    FUN=test_overrepresentation.gene,
                    overrepresentation.frame=overrepresentation.frame)
  overrepresentation.frame$pvalue <- pvalues
  overrepresentation.frame$adj_pvalue <- p.adjust(pvalues, method="BH")
  return(overrepresentation.frame)
}

get_markers.seed <- function(seed, occurrences.population) {
  #' Get genes significantly over-represented in a seed, w.r.t. the population w/ equal sampling effort.
  #'
  #' @param seed: a list with four keys: 'consensus', 'cells', 'clusters' and 'genes'.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #'
  #' @return a vector of characters.
  #'
  overrepresentation.frame <- test_overrepresentation.seed(seed, occurrences.population)
  is_significant <- overrepresentation.frame$adj_pvalue < 0.05
  markers <- rownames(overrepresentation.frame[is_significant, ])
  return(markers)
}

get_markers <- function(seeds, occurrences.population) {
  #' Identify markers w.r.t. their respective seed.
  #'
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'genes'.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #'
  #' @return a nested list with two keys: 'seed' and 'all'.
  #'
  respective_markers <- lapply(X=seeds,
                               FUN=get_markers.seed,
                               occurrences.population=occurrences.population)
  all_markers <- unlist(respective_markers)
  markers <- list(seed=respective_markers, all=all_markers)
  return(markers)
}

add_specific_markers <- function(seeds, markers) {
  #' Add the markers w.r.t. the seeds.
  #'
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param markers: a nested list with two keys: 'seed' and 'all'.
  #'
  #' @return a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #'
  all_markers <- markers[["all"]]
  unspecific_markers <- unique(all_markers[duplicated(all_markers)])
  for (i in 1:length(seeds)) {
    specific_markers <- setdiff(markers$seed[[i]], unspecific_markers)
    seeds[[i]]$specific_markers <- specific_markers
  }
  return(seeds)
}

get_markers.plot <- function(markers, population, params) {
  #' Get an upset plot corresponding to the intersection of marker genes between sets.
  #'
  #' @param markers: a nested list with two keys: 'seed' and 'all'.
  #' @param population: a character.
  #' @param params: a list of parameters, with 'n_HVGs'.
  #'
  #' @return an upset plot.
  #'
  labels <- paste(population, 1:length(markers$seed), sep=".")
  markers.plot <- ggVennDiagram(markers$seed,
                                category.names=labels,
                                force_upset=TRUE,
                                order.set.by="name",
                                order.intersect.by="none",
                                relative_height=1.6,
                                relative_width=0.4)

  ### use the color palette of seeds ###
  ######################################

  # intersection dots
  intersect_colors.dots <- rep(NA, length(markers.plot[[1]]$data$name))
  intersect_colors.dots[1:length(labels)] <- labels
  markers.plot[[1]]$layers[[1]]$aes_params$colour <- NULL
  markers.plot[[1]]$layers[[2]]$aes_params$colour <- NULL
  markers.plot[[1]] <- markers.plot[[1]] +
    aes(colour=intersect_colors.dots) +
    theme(legend.position="none")

  # intersection bars
  intersect_colors.bars <- rep(NA, length(markers.plot[[2]]$data$name))
  intersect_colors.bars[1:length(labels)] <- labels
  markers.plot[[2]] <- markers.plot[[2]] +
    aes(fill=intersect_colors.bars) +
    theme_classic() +
    theme(legend.position="none", axis.text.y=element_text(vjust=0.25),
          panel.grid.major.y=element_line(linewidth=0.5),
          axis.line=element_blank(),
          axis.ticks.x=element_line(colour="#00000000"),
          axis.text.x=element_text(colour="#00000000")) +
    scale_y_continuous(expand=expansion(mult=c(0, .05))) +
    ggtitle("Over-represented HVGs")

  # size bars
  markers.plot[[3]] <- markers.plot[[3]] + aes(fill=name) +
    theme(legend.position="none") +
    geom_vline(xintercept=sqrt(params$n_HVGs), linetype="dashed")

  return(markers.plot)
}

draw_genes <- function(data.loop, seeds, population, params) {
  #' Draw two plots corresponding to:
  #' - a boxplot, with populations as x-axis and number of HVGs measured as y-axis.
  #' - an upsetplot, where each bar corresponds to a set of marker genes.
  #' The marker genes specific to a population are colorized.
  #'
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param seeds: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #' @param population: a character.
  #' @param params: a list of parameters, with 'n_HVGs'.
  #'
  markers.loop <- get_markers(seeds, data.loop$occurrences.loop)
  markers.plot <- get_markers.plot(markers.loop, population, params)

  pdf(file = glue("{params$figures_dir}/{population}_genes.pdf"))
  print(markers.plot)
  dev.off()
}

get_genes <- function(data.loop, seeds, params, population, figures) {
  #' Identify relevant genes to characterize the seeds found.
  #' If too little genes are identified, the iteration is considered seedless.
  #'
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param seeds: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #' @param params: a list of parameters, with 'n_HVGs'.
  #' @param population: a character.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the genes identification.
  #'
  #' @return a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #'

  # get seed-specific markers
  ###########################
  seeds <- add_occurrences_to_seeds(data.loop$ranked_genes.loop, seeds)
  markers.loop <- get_markers(seeds, data.loop$occurrences.loop)
  for (i in 1:length(seeds)) {seeds[[i]]$markers <- markers.loop$seed[[i]]}
  seeds <- add_specific_markers(seeds, markers.loop)
  if (figures) {draw_genes(data.loop, seeds, population, params)}

  # if any seed is poorly characterized, the iteration is fruitless
  #################################################################
  is_informative <- function(seed) {length(seed$markers) >= sqrt(params$n_HVGs)}
  informative_seeds <- Filter(f=is_informative, x=seeds)
  if (length(seeds) != length(informative_seeds)) {
    seeds <- list()
  }
  return(seeds)
}
