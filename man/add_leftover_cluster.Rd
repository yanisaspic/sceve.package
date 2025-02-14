% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_characterized_clusters.R
\name{add_leftover_cluster}
\alias{add_leftover_cluster}
\title{Add a leftover cluster to the list of robust clusters.}
\usage{
add_leftover_cluster(population, robust_clusters, data.iteration, params)
}
\arguments{
\item{population}{a character. It corresponds to the cell population that scEVE will attempt to cluster.}

\item{robust_clusters}{list where every element is a robust pool of cells.
The elements are named lists, with six names:
\code{base_clusters}, \code{cells}, \code{clustering_methods}, \code{label}, \code{markers} and \code{robustness}.}

\item{data.iteration}{a named list, with two names: \code{expression} and \code{SeuratObject}.
They correspond to the scRNA-seq expression matrix of a specific cell population and its SeuratObject.}

\item{params}{a list of parameters (cf. \code{sceve::get_default_parameters()}).}
}
\value{
a list where every element is a pool of cells.
The elements are named lists, with six names:
\code{base_clusters}, \code{cells}, \code{clustering_methods}, \code{label}, \code{markers}, and \code{robustness}.
}
\description{
It corresponds to a group of cells unassigned to any robust cluster.
}
