% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_base_clusters.R
\name{use_Seurat}
\alias{use_Seurat}
\title{Predict clusters with the Seurat method.}
\usage{
use_Seurat(SeuratObject.population, params)
}
\arguments{
\item{SeuratObject.population}{a SeuratObject, on which the function ScaleData()
of Seurat has been applied already.}

\item{params}{a list of parameters (cf. \code{sceve::get_default_parameters()}).}
}
\value{
a named factor that associates each cell to their cluster prediction.
}
\description{
Predict clusters with the Seurat method.
}
