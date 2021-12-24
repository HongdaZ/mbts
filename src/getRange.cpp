#include "getRange.h"

vector<int> getRange( const int *ptr_seg, const int &label,
                      const int *ptr_aidx, 
                      const int &len,
                      const int &nr, const int &nc, const int &ns ) {
  vector<int> c_xyz( 3, 0 );
  // min_x, max_x, min_y, max_y, min_z, max_z
  vector<int> range{ nr + 1, 0, nc + 1, 0, ns + 1, 0 }; 
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] == label ) {
      for( int j = 0; j < 3; ++ j ) {
        c_xyz[ j ] = ptr_aidx[ 3 * i + j ];
        if( range[ 2 * j ] > c_xyz[ j ] ) {
          range[ 2 * j ] = c_xyz[ j ];
        }
        if( range[ 2 * j + 1 ] < c_xyz[ j ] ) {
          range[ 2 * j + 1 ] = c_xyz[ j ];
        }
      }
    }
  }
  return range;
}