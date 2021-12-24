shrink <- function( images ) {
  for( m in 1 : length( images ) ){
    image_shrink <- array(dim = dim( images[[ m ]] ) / 2 )
    idx <- vector( "integer", length = length( image_shrink ) )
    l <- 1
    for( k in 0 : ( dim( image_shrink )[ 3 ] - 1 ) )  {
      for( j in 0 : ( dim( image_shrink )[ 2 ] - 1 ) ) {
        for( i in 0 : ( dim( image_shrink )[ 1 ] - 1 ) ) {
          idx[ l ] <- 240 * 240 * k * 2 + 240 * j * 2 + i * 2 + 1
          l <- l + 1
        }
      }
    }
    image_shrink <- images[[ m ]][ idx ]
    dim( image_shrink ) <- dim( images[[ m ]] ) / 2
    images[[ m ]] <- image_shrink
  }
  images
}
