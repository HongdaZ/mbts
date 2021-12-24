#include "smallRegion.h"

// Find ptr_seg2 voxels within the range
void smallRegion( const vector<int> &range, 
                  const int *ptr_seg2, const int &label2, 
                  const int *ptr_aidx,
                  int *ptr_seg2_copy, const int &len ) {
  vector<int> c_xyz( 3, 0 );
  bool have = true;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg2[ 2 * i ] == label2 ) {
      have = true;
      for( int j = 0; j < 3; ++ j ) {
        c_xyz[ j ] = ptr_aidx[ 3 * i + j ];
        if( c_xyz[ j ] < range[ 2 * j ] ||
            c_xyz[ j ] > range[ 2 * j + 1 ] ) {
          have = false;
        }
      }
      if( have ) {
        ptr_seg2_copy[ 2 * i ] = label2;
      }
    }
  }
}