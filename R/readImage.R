# read the zipped images
readImage <- function( patient ) {
  flair_file <- patient[ 1 ]
  t1ce_file <- gsub( "_flair.nii.gz", "_t1ce.nii.gz", flair_file )
  t2_file <- gsub( "_flair.nii.gz", "_t2.nii.gz", flair_file )
  
  t1ce <- readNIfTI( t1ce_file, reorient = FALSE )@.Data
  flair <- readNIfTI( flair_file, reorient = FALSE )@.Data
  t2 <- readNIfTI( t2_file, reorient = FALSE )@.Data
  
  t1ce[ t1ce == 0 ] <- NaN
  flair[ flair == 0 ] <- NaN
  t2[ t2 == 0 ] <- NaN
  
  res <- list( t1ce = t1ce, flair = flair, t2 = t2 )
  return( res )
}