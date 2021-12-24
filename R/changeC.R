# lbl_flair <- lbl
# lbl_flair[ flair < q_flair[ 4 ] ] <- NA_integer_
changeC <- function( old, image, q ) {
  .Call( "changeC", old, image, q )
}