% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iSpatial.R
\name{infer_v2}
\alias{infer_v2}
\title{Main function of iSpatial}
\usage{
infer_v2(
  spRNA,
  scRNA,
  dims = 1:30,
  k.neighbor = 30,
  infered.assay = "enhanced",
  weighted.KNN = TRUE,
  RNA.weight = 0.5,
  n.core = 8,
  correct.spRNA = FALSE,
  correct.scRNA = FALSE,
  correct.neighbor = 5
)
}
\arguments{
\item{spRNA}{seurat object of spatial transcriptome data.
Should contain normalized data. Run Seurat::NormalizeData or others.}

\item{scRNA}{seurat object of single cell RNA-seq data.
Should contain normalized data. Run Seurat::NormalizeData or others.}

\item{dims}{which dimensions to use when find the nearest neighbors}

\item{k.neighbor}{number of neighbors to use when infer the expression}

\item{infered.assay}{names of output assay in seurat object}

\item{weighted.KNN}{we use a dynamics weight to assign scRNA values to spRNA.
The dynamics weight based on the correlation between scRNA cell and spRNA
cell.}

\item{RNA.weight}{the weight of scRNA-seq expression when genes detected both
in spRNA and scRNA. RNA.weight should be 0 to 1.
The inferred expression = (1 - RNA.weight) * spRNA + RNA.weight * scRNA.}

\item{n.core}{number of CPU cores used to parallel.}

\item{correct.spRNA}{Whether to stabilize expression in spRNA.}

\item{correct.scRNA}{Whether to stabilize expression in scRNA.}

\item{correct.neighbor}{number of nearest neighbors used to correct expr.}
}
\value{
returns a seurat object with inferred expression in infered.assay.
}
\description{
This function integrates scRNA-seq and spatial transcriptome. With limited
genes in merFISH/SeqFISH, iSpatial infer/enhance all genes from scRNA-seq.
Finally, this function generates spatially expression of all genes.
}
