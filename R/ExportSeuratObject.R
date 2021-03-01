#' @title Export a Seurat object in BCS format
#' @param object a Seurat object
#' @param bcs.path path to BCS file
#' @param unique.limit ignore a metadata if it has number of unique labels larger than this number. Default is 100.
#' @param clustering.name name of the metadata for clustering result. Default is "seurat_clusters".
#' @param compression.level an integer ranging from 1 to 9. Higher level creates smaller files but takes more time to create and load. Default is 1.
#' @param author email of the creator. Default is rBCS.
#' @param overwrite if TRUE, overwrite existing file at bcs.path. Default is FALSE.
#' @import Matrix
#' @import rhdf5
#' @importFrom jsonlite write_json
#' @importFrom zip zip
#' @importFrom methods as
#' @export
ExportSeuratObject <- function(
  object,
  bcs.path,
  unique.limit = 100,
  clustering.name = "seurat_clusters",
  compression.level = 1,
  author = "rBCS",
  overwrite = FALSE
) {
  GetSparseMatrix <- function(x) {
    if (class(x)[1] == "dgCMatrix") {
      return(x)
    }
    return(as(x, "sparseMatrix"))
  }

  IsNullMatrix <- function(x) {
    return(is.null(x) || nrow(x) == 0 || ncol(x) == 0)
  }

  # assays: list of Assay class in Seurat
  # keys: a vector of keys sort by priority (descending)
  GetData <- function(assays, assay.name, keys) {
    data <- attr(assays[[assay.name]], keys[1])
    if (IsNullMatrix(data) && length(keys) > 1) {
      warning(paste("Assay", assay.name, "has no", keys[1]))
      data <- GetData(assays, assay.name, keys[-1])
    }
    if (!IsNullMatrix(data)) {
      data <- GetSparseMatrix(data)
    }
    return(data)
  }

  ValidateInput <- function() {
    if (length(object@reductions) == 0) {
      stop("Seurat object has no dimensionality reduction")
    }
    if (!"RNA" %in% names(object@assays)) {
      stop("Seurat object must have an RNA assay")
    }
    if (!all(names(object@assays) %in% c("RNA", "ADT"))) {
      warning("This object has several assays. rBCS only use data from either RNA or ADT.")
    }
    if (file.exists(bcs.path)) {
      if (overwrite) {
        warning(paste(basename(bcs.path), "will be replaced"))
        file.remove(bcs.path)
      } else {
        stop(paste(basename(bcs.path), "already exists. Please use overwrite=TRUE to force replacing."))
      }
    }
  }

  GetExpressionData <- function(object) {
    omics <- names(object@assays)
    counts <- GetData(object@assays, "RNA", c("counts", "data"))
    norms <- GetData(object@assays, "RNA", "data")
    feature.type <- rep("RNA", nrow(norms))
    feature.name <- rownames(norms)
    if ("ADT" %in% omics) {
      counts <- Matrix::rbind2(counts, GetData(object@assays, "ADT", c("counts", "data")))
      norms <- Matrix::rbind2(norms, GetData(object@assays, "ADT", "data"))
      feature.type <- c(feature.type, rep("ADT", nrow(norms) - length(feature.type)))
      feature.name <- rownames(norms)
      feature.name[feature.type == "ADT"] <- paste("ADT", feature.name[feature.type == "ADT"], sep="-")
    }
    return(list(counts=counts, norms=norms, feature_type=feature.type))
  }

  GetDimredData <- function(object) {
    data <- lapply(object@reductions, function(dimred) {
      info <- list(
        param = list(omics=dimred@assay.used),
        coords = as.matrix(dimred@cell.embeddings)
      )
      if (ncol(info$coords) > 3) {
        info$coords <- info$coords[, 1:3]
      }
      return(info)
    })
    names(data) <- names(object@reductions)
    return(data)
  }

  ValidateInput()

  Meow("Initializing...")
  hash <- uuid::UUIDgenerate()
  old.wd <- getwd()
  tmp.wd <- dirname(bcs.path)
  dir.create(tmp.wd, recursive=TRUE, showWarnings=FALSE)
  setwd(tmp.wd) # easier for zipping
  on.exit(setwd(old.wd))
  dir.create(file.path(hash, "main"), recursive=TRUE)

  Meow("Extracting expressions...")
  expr.data <- GetExpressionData(object)
  dimred.data <- GetDimredData(object)
  WriteStudy(expr.data, object@meta.data, dimred.data, hash,
      clustering.name, author, unique.limit)

  Meow("Compressing data...")
  zip::zip(bcs.path, hash, compression_level=compression.level)
  unlink(hash, recursive=TRUE, force=TRUE)
  return(TRUE)
}

Meow <- function(...) {
  prefix <- paste0("[", Sys.time(), "]")
  cat(prefix, ..., "\n")
}
