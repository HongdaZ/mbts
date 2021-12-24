# Check is points are inside the convex hull
inSide <- function( p, poly ) {
  .Call( "inSide", p, poly )
}