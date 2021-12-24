#include "findOutLabel.h"

// find an outlier label for current voxel
void findOutLabel( int &out_label, const int &idx, const int *ptr_seg, 
                   const vector<int> &outl_labels ) {
  // current label is not outlier
  if( ptr_seg[ 2 * ( idx - 1 ) ] < 1 ) {
    int i = 0;
    while( i < outl_labels.size() && outl_labels[ i ] != 0 ) {
      ++ i;
    }
    out_label = i + 1;
  } else {
    out_label = ptr_seg[ 2 * ( idx - 1 ) ];
  }
  return;
}