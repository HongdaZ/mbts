#' A Three-stage Procedure for Multimodal Brain Tumor Segmentation
#' 
#' \code{mbts} analyzes T1ce, FLAIR and T2 images and segments the images 
#' into brain tissues.
#' @author Hongda Zhang
#' @useDynLib mbts
#' @importFrom SimpleITK ReadImage Cast N4BiasFieldCorrectionImageFilter WriteImage
#' @param infile Full path of original MR image.
#' @param outdir The directory for output file.
#' @param iter Iterations to run.
#' @param res Level of resolution.
#' @export 
biasFieldCorrect <-  function(infile, outdir, iter = 100, res = 3) {
  ## checks
  stopifnot(file.exists(infile))
  stopifnot(dir.exists(outdir))
  
  ## write result into file in 'outdir' with same base name as 'infile'
  outfile <- file.path(outdir, basename(infile))
  
  ## Read data
  input <- ReadImage(infile)
  input <- Cast(input, "sitkFloat32")
  
  ## Bias field correction
  corrector <- N4BiasFieldCorrectionImageFilter()
  corrector$SetMaximumNumberOfIterations(rep(iter, res))
  output <- corrector$Execute(input)
  
  ## Write the image into the output folder
  WriteImage(output, outfile)
}