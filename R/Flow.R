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

samplingCells <- function(fs, numCells, isEven = T) {
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
	for (i in 1:n) {
		fcs <- fs[[i]]
		fs[[i]] <- fcs[sample(1:nrow(fcs), cellsPerSample[i]), ]
	}
	return(fs)
}
