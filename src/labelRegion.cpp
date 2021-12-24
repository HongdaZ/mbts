#include "labelRegion.h"


// the voxel in region is labelled as 3 in ptr_seg[ 2, ]
void labelRegion( const vector<int> &region,
                  int *ptr_seg ) {
  
  for( int i = 0; i < region.size(); ++ i ) {
    ptr_seg[ 2 * region[ i ] - 1 ] = 3;
  }
  return;
}
// change ptr_seg[ 2, region ] back
void recoverLabel( const vector<int> &region,
                   int *ptr_seg ) {
  
  for( int i = 0; i < region.size(); ++ i ) {
    ptr_seg[ 2 * region[ i ] - 1 ] = 0;
  }
  return;
}