#' @title Write a BBrowser study
#' @param expr.data list of counts, norms, and feature_type
#' @param metadata a data.frame of metadata
#' @param dimred.data a named list of dimensionality reductions
#' @param study.path path to study directory
#' @param clustering.name name of metadata that indicates clustering result
#' @param author email of the creator
#' @importFrom jsonlite write_json
#' @importFrom uuid UUIDgenerate
WriteStudy <- function(
  expr.data, metadata, dimred.data, study.path, clustering.name, author, unique.limit
) {
  CreateCommit <- function() {
    return(list(
      hash_id = uuid::UUIDgenerate(),
      created_at = as.numeric(Sys.time()) * 1000,
      created_by = author,
      description = "Converted with rBCS"
    ))
  }

  WriteH5Matrix <- function(matrix, h5, group, write.labels=TRUE) {
    rhdf5::h5createGroup(h5, group)
    rhdf5::h5write(matrix@x, h5, paste0(group, "/data"))
    rhdf5::h5write(matrix@p, h5, paste0(group, "/indptr"))
    rhdf5::h5write(matrix@i, h5, paste0(group, "/indices"))
    rhdf5::h5write(dim(matrix), h5, paste0(group, "/shape"))
    if (write.labels) {
      rhdf5::h5write(colnames(matrix), h5, paste0(group, "/barcodes"))
      rhdf5::h5write(rownames(matrix), h5, paste0(group, "/features"))
    }
  }

  WriteExpressionData <- function(expr.data, study.path) {
    Meow("Writing column sums...")
    h5.path <- file.path(study.path, "main", "matrix.hdf5")
    rhdf5::h5createFile(h5.path)
    h5 <- rhdf5::H5Fopen(h5.path)
    rhdf5::h5createGroup(h5, "colsum")
    raw.sum <- Matrix::colSums(expr.data$counts)
    rhdf5::h5write(raw.sum, h5, "colsum/raw")
    rhdf5::h5write(log2(raw.sum + 1), h5, "colsum/log")
    rhdf5::h5write(Matrix::colSums(expr.data$norms), h5, "colsum/lognorm")
    Meow("Writing matrices...")
    WriteH5Matrix(expr.data$counts, h5, "bioturing")
    WriteH5Matrix(Matrix::t(expr.data$counts), h5, "countsT")
    WriteH5Matrix(Matrix::t(expr.data$norms), h5, "normalizedT")
    rhdf5::H5Fclose(h5)
    writeLines(rownames(expr.data$norms), file.path(study.path, "main", "genes.tsv"))
    writeLines(colnames(expr.data$norms), file.path(study.path, "main", "barcodes.tsv"))
  }

  WriteMetadata <- function(metadata, study.path, clustering.name) {
    Meow("Writing metadata...")
    meta.dir <- file.path(study.path, "main", "metadata")
    dir.create(meta.dir)
    meta.list <- list(version=1)
    meta.list$content <- lapply(seq(ncol(metadata)), function(i) {
      info <- list(
        id = uuid::UUIDgenerate(),
        name = colnames(metadata)[i],
        type = if (is.numeric(metadata[[i]])) "numeric" else "category",
        clusters = as.character(metadata[[i]]),
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
        info$clusters[is.na(info$clusters)] <- "Unassigned"
        info$clusterName <- c("Unassigned", unique(info$clusters))
        if (length(info$clusterName) > unique.limit) {
          warning(paste("Bad metadata -", info$name))
          info$bad <- TRUE
          return(info)
        }
        info$clusterLength <- sapply(info$clusterName, function(x) sum(info$clusters == x))
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
  }

  WriteDimred <- function(dimred.data, study.path) {
    Meow("Writing cell embeddings...")
    dimred.dir <- file.path(study.path, "main", "dimred")
    dir.create(dimred.dir)
    dimred <- list(version=1)
    dimred$data <- lapply(seq(length(dimred.data)), function(i) {
      info <- list(
        name = names(dimred.data)[i],
        id = uuid::UUIDgenerate(),
        param = dimred.data[[i]]$param,
        coords = dimred.data[[i]]$coords,
        history = list(CreateCommit())
      )
      info$size <- dim(info$coords)
      if (dimred.data[[i]]$type == "dimred") {
        jsonlite::write_json(info, file.path(dimred.dir, info$id), auto_unbox=TRUE)
      }
      info$coords <- NULL # meta does not need this info
      return(info)
    })
    dimred$data <- dimred$data[sapply(dimred.data, function(x) x$type == "dimred")]
    names(dimred$data) <- sapply(dimred$data, function(x) x$id)
    dimred$default <- tail(names(dimred$data), 1)
    jsonlite::write_json(dimred, file.path(dimred.dir, "meta"), auto_unbox=TRUE)
  }

  WriteLowDimred <- function(dimred.data, study.path) {
    Meow("Writing intermediate embeddings...")
    h5.path <- file.path(study.path, "main", "pca_result.hdf5")
    h5 <- rhdf5::H5Fcreate(h5.path)
    for (i in seq_along(dimred.data)) {
      coords <- dimred.data[[i]]$coords
      if (dimred.data[[i]]$type == "low_dimred") {
        if (class(coords)[1] == "matrix") {
          rhdf5::h5write(coords, h5.path, names(dimred.data)[i])
        } else if (class(coords)[1] == "dgCMatrix") {
          WriteH5Matrix(coords, h5.path, names(dimred.data)[i], write.labels=FALSE)
        }
      }
    }
    rhdf5::H5Fclose(h5)
  }

  WriteRunInfo <- function(study.path, omics, n.cell, title) {
    Meow("Writing general information...")
    run.info <- list(
      species = "human",
      omics = as.list(omics), # Must be list
      hash_id = basename(study.path),
      version = 16,
      n_cell = n.cell,
      addon = "SingleCell",
      n_batch = 1,
      platform = "unknown",
      title = "Untitled study",
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
    jsonlite::write_json(run.info, file.path(study.path, "run_info.json"), auto_unbox=TRUE)
  }

  dir.create(file.path(study.path, "main"), recursive=TRUE)
  WriteExpressionData(expr.data, study.path)
  WriteMetadata(metadata, study.path, clustering.name)
  WriteDimred(dimred.data, study.path)
  WriteLowDimred(dimred.data, study.path)
  WriteRunInfo(study.path, unique(expr.data$feature.type), ncol(expr.data$norms))
}
