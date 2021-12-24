# img[ label != 0 ] <- NaN
changeA <- function( img, label ) {
  .Call( "changeA", img, label )
}