#include "excldVoxel.h"

// Exclude voxels with ptr_seg[ 2 * i ] == label (region[ i ] = 0)
void excldVoxel( vector<int> &region, const int *ptr_seg, 
                 const int &label ) {
  int curr_idx;
  for( vector<int>::iterator it = region.begin(); it != region.end();
  ++ it ) {
    curr_idx = *it;
    if( curr_idx != 0 ) {
      if( ptr_seg[ 2 * ( curr_idx - 1 ) ] == label ) {
        *it = 0;
      }
    }
  }
}