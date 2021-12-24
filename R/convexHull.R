# Find the convex hull
convexHull <- function( points ) {
  .Call( "convexHull", points )
}