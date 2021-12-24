# Return voxel index vector, neighbor index matrix, intensity vector
# and neighbor intensity matrix
indexMat <- function(  img, label ) {
  .Call( "indexMat", img, label )
}