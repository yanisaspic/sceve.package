% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_characterized_clusters.R
\name{draw_marker_genes}
\alias{draw_marker_genes}
\title{Get an upset-plot representing the marker genes predicted in each meta-cluster.}
\usage{
draw_marker_genes(characterized_clusters)
}
\arguments{
\item{characterized_clusters}{a list where every element is a characterized pool of cells.
The elements are named lists, with six names:
\code{base_clusters}, \code{cells}, \code{clustering_methods}, \code{label}, \code{markers} and \code{robustness}.}
}
\value{
a plot.
}
\description{
Get an upset-plot representing the marker genes predicted in each meta-cluster.
}
