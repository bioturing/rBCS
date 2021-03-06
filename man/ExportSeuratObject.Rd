% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportSeuratObject.R
\name{ExportSeuratObject}
\alias{ExportSeuratObject}
\title{Export a Seurat object in BCS format}
\usage{
ExportSeuratObject(
  object,
  bcs.path,
  unique.limit = 100,
  clustering.name = "seurat_clusters",
  compression.level = 1,
  author = "rBCS",
  overwrite = FALSE
)
}
\arguments{
\item{object}{a Seurat object}

\item{bcs.path}{path to BCS file}

\item{unique.limit}{ignore a metadata if it has number of unique labels larger than this number. Default is 100.}

\item{clustering.name}{name of the metadata for clustering result. Default is "seurat_clusters".}

\item{compression.level}{an integer ranging from 1 to 9. Higher level creates smaller files but takes more time to create and load. Default is 1.}

\item{author}{email of the creator. Default is rBCS.}

\item{overwrite}{if TRUE, overwrite existing file at bcs.path. Default is FALSE.}
}
\description{
Export a Seurat object in BCS format
}
