## Initialize data for post processing
initPost <- function( t1ce_image, flair_image, t2_image,
                      t1ce_intst ) {
  ## Return vector index, neighbor index, array index, t1ce_seg,
  ## flair_seg, t2_seg
  .Call( "initPost", t1ce_image, flair_image, t2_image,
         t1ce_intst )
}