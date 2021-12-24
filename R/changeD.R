# label[ label == 1 | label == 4 ] <- NA_integer_
changeD <- function( label, a, b ) {
  .Call( "changeD", label, a, b )
}