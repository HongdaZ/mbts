getB <- function( m, a ) {
  ( a + 1 ) * ( m[ length( m ) ] - max( m[ - length( m ) ] ) ) * 3.5
}