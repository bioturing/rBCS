readFlowSet <- function(path) {
	### Read each fcs file into a dataframe
	require(flowCore)
	fcs.files <- list.files(path)
	fcs.files <- fcs.files[grep(".fcs$", fcs.files)]
	fcs.files <- lapply(fcs.files, function(f) {
		return(file.path(path, f))
	})
	sample.names <- lapply(fcs.files, basename)
	fs <- lapply(fcs.files, read.FCS)
	names(fs) <- sample.names
	message(paste("Read", length(fcs.files), "samples"))
	print(as.character(fcs.files))
	return(fs)
}

getMetadata <- function(fs) {
	expName <- as.character(lapply(fs, function(fcs) {
						return(fcs@description$`EXPERIMENT NAME`)
				}))
	plateName <- as.character(lapply(fs, function(fcs) {
						return(fcs@description$`PLATE NAME`)
				}))
	channelDetails <- lapply(fs, function(fcs) {
						return(fcs@parameters@data)
				})
	### The most important one
	exprsFrame <- lapply(fs, function(fcs) {
						return(data.frame(fcs@exprs))
				})
	message("Extract the expName, plateName, channelDetails, and expression table")
	return(list(expName=expName, 
				plateName=plateName, 
				channelDetails=channelDetails,
				exprsFrame=exprsFrame))
}

renameFCS <- function(fs, channelDetails) {
	### convert space to underscore and handle <NA>
	modifyName <- function(i, j) {
		fcs <- channelDetails[[i]]
		### i: index of fsc, j: index of feature in each fcs
		if (is.na(fcs$desc[j])) {
			return(gsub(" |-", "_", fcs$name[j]))
		}
		return(gsub(" |-", "_", paste0(fcs$desc[j], "_(",fcs$name[j], ")")))
	}
	xcolnames <- lapply(1:length(fs), function(i) unlist(lapply(1:nrow(channelDetails[[i]]), function(j) {
				return(modifyName(i, j))
			}
	)))
	message("Rename channels...")

	for (i in 1:length(fs)) {
		colnames(fs[[i]]) <- as.character(xcolnames[[i]])
	}
	return(fs)
}

getCellsPerSample <- function(fs, numCells, isEven = T) {
	n <- length(fs)
	if (isEven) {
		message("Even sampling...")
		targetedNumber <- min(floor(numCells / n), as.integer(lapply(fs, nrow)))
		message(paste(targetedNumber, "is the targeted number of cells for each sample"))
		cellsPerSample <- rep(targetedNumber, n)
		print(cellsPerSample)
	} else {
		message("Non-even sampling...")
		nCellsPerSample <- as.integer(lapply(fs, nrow))
		origTotalCells <- sum(nCellsPerSample)
		targetedRatios <- numCells / origTotalCells
		cellsPerSample <- floor(targetedRatios * nCellsPerSample)
		print(cellsPerSample)
	}
	return(cellsPerSample)
}

samplingCells <- function(fs, cellsPerSample) {
	for (i in 1:length(fs)) {
		fcs <- fs[[i]]
		fs[[i]] <- fcs[sample(1:nrow(fcs), cellsPerSample[i]), ]
	}
	return(fs)
}

mergeSample <- function(fs) {
	sharedChannels <- Reduce(intersect, lapply(fs, colnames))
	message(paste("There is total", length(sharedChannels), "shared channels in the data"))
	if (length(sharedChannels) == 0) {
		stop("Cannot merge all samples together because nothing is shared among them!")
	}
	print(sharedChannels)
	fs <- lapply(fs, function(fcs) fcs[, sharedChannels])
	message("Merging frames...")
	fs <- Reduce(rbind, fs)
	message("All frames merged")
	return(fs)
}

removeScatter <- function(fss) {
	isScatterChannels <- grepl("SSC", colnames(fss)) | grepl("FSC", colnames(fss))
	fss <- fss[, !isScatterChannels]
	return(fss)
}

makeNonNegative <- function(x) {
	if (min(x) < 0) {
		x <- x + (0 - min(x))
	}
}

addMetadata <- function(fssMeta, cellsPerSample, obj) {
	n <- length(cellsPerSample)
	if (!is.null(fssMeta$expName)) {
		obj$Experiment_Name <- as.character(sapply(1:n, function(i) {
							rep(fssMeta$expName[i], cellsPerSample[i])
				}))
	}
	if (!is.null(fssMeta$plateName)) {
		obj$Plate_Name <- as.character(sapply(1:n, function(i) {
							rep(fssMeta$plateName[i], cellsPerSample[i])
				}))
	}
	obj$Sample_Name <- as.character(sapply(1:n, function(i) {
							rep(fssMeta$sampleName[i], cellsPerSample[i])
				}))
	return(obj)
}

DoingSeuratStuff <- function(obj) {
	require(Seurat)
	obj <- ScaleData(obj)
	obj <- RunPCA(obj, features = rownames(obj))
	obj@reductions$xpca <- obj@reductions$pca
	obj@reductions$pca <- NULL
	return(obj)
}

ExportFCS <- function(FCSpath,
					  BCSpath,
					  numCells = 1000000, 
					  evenSampling = T, 
					  removingScatter = T,
					  shiftNegative = T) {
	require(Seurat)
	fs <- readFlowSet(FCSpath)
	m <- getMetadata(fs)
	fs <- m$exprsFrame
	m$exprsFrame <- NULL
	m$sampleName <- names(fs)
	fs <- renameFCS(fs, m$channelDetails)
	cellsPerSample <- getCellsPerSample(fs, numCells, evenSampling)
	fs <- samplingCells(fs, cellsPerSample)
	fss <- mergeSample(fs)
	if (removingScatter) {
		fss <- removeScatter(fss)
	}
	cellTime <- fss$Time
	fss$Time <- NULL
	if (removingScatter) {
		message("Removed Time and scattering channels")
	} else {
		message("Removed Time channel")
	}
	if (shiftNegative) {
		fss <- apply(fss, 2, makeNonNegative)
		message("Each channel was shifted right by the min of its values")
	} else {
		message("Negative values were converted to zeroes")
		fss <- apply(fss, 2, function(x) {
					 x[x < 0] <- 0
					 return(x)
		})
	}
	fss <- t(fss)
	obj <- CreateSeuratObject(counts = fss, project = "fcs_integration")
	obj <- addMetadata(m, cellsPerSample, obj)
	message("Converted fss to the Seurat object")
	obj <- DoingSeuratStuff(obj)
	ExportSeuratObject(obj, BCSpath)
}

