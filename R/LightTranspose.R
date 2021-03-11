#' @title Tranpose a CSC sparse matrix with low memory and high speed 
#' @param input_path Input path of the hdf5 matrix 
#' @param output_path Path of the hdf5 matrix 
#' @param chunk_size Divide the matrix into equal chunks. It should be the square root of number of columns
#' @export
LightTranspose <- function(
  input.path,
  output.path,
  chunk_size = 100
  ) {
  
  light.tranpose.bin <- file.path(find.package("rBCS"), "lib", "light_transpose.pyx.dist", 
                                  "light_transpose.pyx")
  if (! file.exists(light.tranpose.bin)) {
    stop("the backend library of light transpose is not available. Your rBCS installation might be broken. Please try to install it again.")
  }
  message("Now transpose the matrix...")
  ret <- system(paste(light.tranpose.bin, input.path, output.path, chunk_size))
  if (ret)
    stop("Some error happened, Please post an issue here to report it: https://github.com/bioturing/rBCS/issues")
  message(paste("Transpose successfully, your output hdf5 file is located at", output.path ))
}
