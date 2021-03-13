# Detect the OS: Linux, Windows, or Darwin
os <- Sys.info()[['sysname']]

# The installation destination. Please use the same way to load the binary
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
if (! dir.exists(dest)) {
	dir.create(dest, recursive = TRUE, showWarnings = TRUE)
}

if (os == "Linux") {
	files <- file.path(R_PACKAGE_SOURCE, "src/bcs_transpose.pyx.dist.linux")
} else if (os == "Windows") {
	files <- file.path(R_PACKAGE_SOURCE, "src/bcs_transpose.pyx.dist.win")
} else if (os == "Darwin") {
	files <- file.path(R_PACKAGE_SOURCE, "src/bcs_transpose.pyx.dist.mac")
}

# Copy file
message("Copying files from ", files, " to ", dest)
message(paste(list.files(files), sep = "\n"))
file.copy(files, dest, overwrite = TRUE, recursive = T, )


