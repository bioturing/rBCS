#' @title Tranpose a CSC sparse matrix with low memory and high speed 
#' @param hdf5.path Input path of the hdf5 matrix. Will expand this matrix file
#' @param chunk_size Divide the matrix into equal chunks. It should be the square root of number of columns
#' @export
LightTranspose <- function(
  hdf5.path,
  chunk_size = 100
  ) {
  
  GetBinary <- function() {
    os <- Sys.info()[['sysname']]
    r_arch <-.Platform $ r_arch
    
    if (os == "Linux") {
      transpose.bin <- system.file("libs", r_arch, "bcs_transpose.pyx.dist.linux/bcs_transpose.pyx", package = "rBCS")
    } else if (os == "Windows") {
      transpose.bin <- system.file("libs", r_arch, "bcs_transpose.pyx.dist.win/bcs_transpose.pyx.exe", package = "rBCS")
      
    } else if (os == "Darwin") {
      transpose.bin <- system.file("libs", r_arch, "bcs_transpose.pyx.dist.mac/bcs_transpose.pyx", package = "rBCS")
    }
    return(transpose.bin)
  }
  
  transpose.bin <- GetBinary()
  
  if (transpose.bin == "") {
    stop("the backend library of the bcs file transpose is not available. Your rBCS installation might be broken. Please try to install it again.")
  }
  message("Now transpose the matrix...")
  cmd <- paste0('"',transpose.bin,'"', ' ' , '"', hdf5.path, '"', ' ', chunk_size)
  message(cmd)
  ret <- system(cmd)
  message(ret)
  if (ret)
    stop("Some error happened, Please post an issue here to report it: https://github.com/bioturing/rBCS/issues")
  message(paste("Transpose successfully, your output hdf5 file is located at", hdf5.path ))
}
