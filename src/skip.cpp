#include "skip.h"

// whether skip current voxel
bool skip( const int &idx, const int *ptr_seg, const int *ptr_nidx ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  if( curr_label != 0 ) {
    return false;
  } else {
    int nbr_idx;
    for( int i = 0; i < 6; ++ i ) {
      nbr_idx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( nbr_idx != NA_INTEGER ) {
        if( ptr_seg[ 2 * ( nbr_idx - 1 ) ] != 0 ) {
          return false;
        }
      }
    }
    return true;
  }
}