% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_data.R
\name{get_expression.selected}
\alias{get_expression.selected}
\title{Get a scRNA-seq dataset of raw count expression, with selected genes.}
\usage{
get_expression.selected(expression, params)
}
\arguments{
\item{expression}{a scRNA-seq dataset of raw count expression, without selected genes.
Its rows are genes and its columns are cells.}

\item{params}{a list of parameters (cf. \code{sceve::get_default_parameters()}).}
}
\value{
a scRNA-seq dataset of raw count expression, with selected genes.
}
\description{
Get a scRNA-seq dataset of raw count expression, with selected genes.
}
