#' @title Export a Seurat object in BCS format
#' @param object a Seurat object
#' @param bcs.path path to BCS file
#' @param unique.limit ignore a metadata if it has number of unique labels larger than this number. Default is 100.
#' @param clustering.name name of the metadata for clustering result. Default is "seurat_clusters".
#' @param compression.level an integer ranging from 1 to 9. Higher level creates smaller files but takes more time to create and load. Default is 1.
#' @param author email of the creator
#' @import Matrix
#' @import rhdf5
#' @importFrom jsonlite write_json
#' @importFrom uuid UUIDgenerate
#' @importFrom zip zip
#' @importFrom methods as
#' @export
ExportSeuratObject <- function(
  object,
  bcs.path,
  unique.limit = 100,
  clustering.name = "seurat_clusters",
  compression.level = 1,
  author = "rBCS"
) {
  GetSparseMatrix <- function(x) {
    if (class(x)[1] == "dgCMatrix") {
      return(x)
    }
    return(as(x, "sparseMatrix"))
  }

  IsNullMatrix <- function(x) {
    return(is.null(x) || (nrow(x) == 0 || ncol(x) == 0))
  }

  # assays: list of Assay class in Seurat
  # key: a vector of key sort by priority (descending)
  GetData <- function(assays, assay.name, key) {
    data <- attr(assays[[assay.name]], key[1])
    if (IsNullMatrix(data) && length(key) > 1) {
      Meow("[WARNING] Assay", assay.name, "has no", key[1])
      data <- GetData(assays, assay.name, key[-1])
    }
    if (!IsNullMatrix(data)) {
      data <- GetSparseMatrix(data)
    }
    return(data)
  }

  WriteH5Matrix <- function(m, group, need.transpose=FALSE) {
    rhdf5::h5createGroup(h5, group)
    if (need.transpose) {
      m <- Matrix::t(m)
    }
    rhdf5::h5write(m@x, h5, paste0(group, "/data"))
    rhdf5::h5write(m@p, h5, paste0(group, "/indptr"))
    rhdf5::h5write(m@i, h5, paste0(group, "/indices"))
    rhdf5::h5write(dim(m), h5, paste0(group, "/shape"))
    rhdf5::h5write(colnames(m), h5, paste0(group, "/barcodes"))
    if (need.transpose) {
      rhdf5::h5write(rownames(m), h5, paste0(group, "/features"))
    } else {
      rhdf5::h5write(feature.name, h5, paste0(group, "/features"))
      rhdf5::h5write(feature.type, h5, paste0(group, "/feature_type"))
    }
  }

  CreateCommit <- function() {
    return(list(
      hash_id = uuid::UUIDgenerate(),
      created_at = as.numeric(Sys.time()) * 1000,
      created_by = author,
      description = "Converted with rBCS"
    ))
  }

  ValidateObject <- function() {
    if (length(object@reductions) == 0) {
      stop("Seurat object has no dimensionality reduction")
    }
    if (!"RNA" %in% names(object@assays)) {
      stop("Seurat object must have an RNA assay")
    }
    if (!all(names(object@assays), c("RNA", "ADT"))) {
      Meow("[WARNING] This object has several assays. rBCS only use data from either RNA or ADT.")
    }
  }

  ValidateObject()

  Meow("Initializing...")
  hash <- uuid::UUIDgenerate()
  old.wd <- getwd()
  tmp.wd <- dirname(bcs.path)
  dir.create(tmp.wd, recursive=TRUE, showWarnings=FALSE)
  setwd(tmp.wd) # easier for zipping
  on.exit(setwd(old.wd))
  dir.create(file.path(hash, "main"), recursive=TRUE)

  Meow("Extracting expressions...")
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

  Meow("Writing column sums...")
  h5.path <- file.path(hash, "main", "matrix.hdf5")
  rhdf5::h5createFile(h5.path)
  h5 <- rhdf5::H5Fopen(h5.path)
  rhdf5::h5createGroup(h5, "colsum")
  raw.sum <- Matrix::colSums(counts)
  rhdf5::h5write(raw.sum, h5, "colsum/raw")
  rhdf5::h5write(log2(raw.sum + 1), h5, "colsum/log")
  rhdf5::h5write(Matrix::colSums(norms), h5, "colsum/lognorm")

  Meow("Writing matrices...")
  WriteH5Matrix(counts, "bioturing")
  WriteH5Matrix(counts, "countsT", need.transpose=TRUE)
  WriteH5Matrix(norms, "normalizedT", need.transpose=TRUE)
  rhdf5::H5Fclose(h5)
  writeLines(feature.name, file.path(hash, "main", "genes.tsv"))
  writeLines(colnames(norms), file.path(hash, "main", "barcodes.tsv"))
  
  Meow("Writing metadata...")
  meta.dir <- file.path(hash, "main", "metadata")
  dir.create(meta.dir)
  meta.list <- list(version=1)
  meta.list$content <- lapply(seq(ncol(object@meta.data)), function(i) {
    info <- list(
      id = uuid::UUIDgenerate(),
      name = colnames(object@meta.data)[i],
      type = if (is.numeric(object@meta.data[[i]])) "numeric" else "category",
      clusters = object@meta.data[[i]],
      clusterName = "NaN",
      clusterLength = 0,
      history = list(CreateCommit())
    )

    # Use existing clustering result
    if (info$name == clustering.name) {
      info$type <- "category"
      info$name <- "Graph-based clusters"
      info$id <- "graph_based"
      if (!any(is.na(as.numeric(info$clusters)))) {
        info$clusters <- paste("Cluster", as.numeric(info$clusters))
      }
    }

    # Category counting
    if (info$type == "category") {
      info$clusterName <- c("Unassigned", as.character(unique(info$clusters)))
      if (length(info$clusterName) > unique.limit) {
        Meow("WARNING: Bad metadata -", info$name)
        info$bad <- TRUE
        return(info)
      }
      info$clusterLength <- c(0, as.numeric(table(info$clusters)))
      info$clusters <- match(as.character(info$clusters), info$clusterName) - 1 # base-0
    }

    jsonlite::write_json(info, file.path(meta.dir, paste0(info$id, ".json")), auto_unbox=TRUE)
    info$clusters <- NULL # metalist.json does not need this info
    return(info)
  })
  meta.list$content <- meta.list$content[sapply(meta.list$content, function(x) is.null(x$bad))]
  names(meta.list$content) <- sapply(meta.list$content, function(x) x$id)
  
  # Add a fake graph-based result. TODO: Recognize graph-based result from object
  if (is.null(meta.list$content$graph_based)) {
    graph <- list(
      id = "graph_based",
      name = "Graph-based clusters",
      type = "category",
      clusters = rep(1, ncol(norms)),
      clusterName = c("Unassigned", "Cluster 1"),
      clusterLength = c(0, ncol(norms)),
      history = list(CreateCommit())
    )
    jsonlite::write_json(graph, file.path(meta.dir, "graph_based.json"), auto_unbox=TRUE)
    graph$clusters <- NULL
    meta.list$content$graph_based <- graph
  }

  jsonlite::write_json(meta.list, file.path(meta.dir, "metalist.json"), auto_unbox=TRUE)
  
  Meow("Writing cell embeddings...")
  dimred.dir <- file.path(hash, "main", "dimred")
  dir.create(dimred.dir)
  dimred <- list(version=1)
  dimred$data <- lapply(seq(length(object@reductions)), function(i) {
    info <- list(
      name = names(object@reductions)[i],
      id = uuid::UUIDgenerate(),
      param = list(omics=object@reductions[[i]]@assay.used),
      coords = as.matrix(object@reductions[[i]]@cell.embeddings),
      history = list(CreateCommit())
    )
    if (ncol(info$coords) > 3) {
      info$coords <- info$coords[, 1:3]
    }
    info$size <- dim(info$coords)
    jsonlite::write_json(info, file.path(dimred.dir, info$id), auto_unbox=TRUE)
    info$coords <- NULL # meta does not need this info
    return(info)
  })
  names(dimred$data) <- sapply(dimred$data, function(x) x$id)
  jsonlite::write_json(dimred, file.path(dimred.dir, "meta"), auto_unbox=TRUE)
  
  Meow("Writing general information...")
  run.info <- list(
    species = "human",
    omics = as.list(unique(feature.type)), # Must be list
    hash_id = hash,
    version = 16,
    n_cell = ncol(object),
    addon = "SingleCell",
    n_batch = 1,
    platform = "unknown",
    title = object@project.name,
    unit = "umi",
    author = list(),
    abstract = "",
    ana_setting = list(inputType="bcs"),
    remote = list(),
    history = list(CreateCommit()),
    tag = list(),
    shareTag = list(),
    modified_date = as.numeric(Sys.time()) * 1000,
    created_date = as.numeric(Sys.time()) * 1000
  )
  jsonlite::write_json(run.info, file.path(hash, "run_info.json"), auto_unbox=TRUE)

  Meow("Compressing data...")
  zip::zip(bcs.path, hash, compression_level=compression.level)
  unlink(hash, recursive=TRUE, force=TRUE)
  return(TRUE)
}

Meow <- function(...) {
  prefix <- paste0("[", Sys.time(), "]")
  cat(prefix, ..., "\n")
}
