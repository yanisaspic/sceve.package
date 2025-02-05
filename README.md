# scEVE: An Alternative Framework to scRNA-seq Ensemble Clustering

Single-cell RNA-sequencing, or scRNA-seq, is a technology developed to measure the individual transcriptomes of cells composing a tissue. In the past decade only, that technology has motivated the development of hundreds of clustering methods that attempt to identify cell populations with similar transcriptomes. Because each method relies on specific hypotheses, their predictions can vary drastically, and relying on a single clustering method is ill-advised.

To address that issue, ensemble clustering algorithms integrate multiple clustering methods and minimize the differences of their predictions to generate a unique result. While that approach is sensible, multiple key challenges remain unadressed in single-cell clustering; namely, ensemble clustering methods have yet to predict clusters at multiple resolutions, and to quantify the uncertainty of their predictions. In this repository, we propose scEVE, an original scRNA-seq ensemble clustering framework, that addresses these two challenges by describing the differences between clustering results, instead of minimizing them.

scEVE is maintained by Asloudj Yanis [yanis.asloudj@u-bordeaux.fr].

## Installation

You can install scEVE from Github with:

```{r}
install.packages("devtools")
devtools::install_github("yanisaspic/sceve.package")
```

## Overview of the scEVE framework

A complete overview of the framework is available in its associated vignette.
