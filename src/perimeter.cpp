#include "Rinternals.h"


#include <vector>

#include "perimeter.h"

using std::vector;

int perimeter( const int *ptr_seg, const int &label, 
               const int &len, const int *ptr_nidx,
               const int &plane ) {
  int count = 0, nidx = 0;
  const vector<int> plane_idx{ 0, 1, 4, 5,
                               0, 2, 3, 5, 
                               1, 2, 3, 4 };
  const vector<int> curr_plane{ plane_idx[ 4 * plane ],
                                plane_idx[ 4 * plane + 1 ],
                                plane_idx[ 4 * plane + 2 ],
                                plane_idx[ 4 * plane + 3 ] };
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] != label ) {
      for( int j = 0; j < 4; ++ j ) {
        nidx = ptr_nidx[ 6 * i + curr_plane[ j ] ];
        if( nidx == NA_INTEGER ) {
          continue;
        } else {
          if( ptr_seg[ 2 * ( nidx - 1 ) ] == label ) {
            ++ count;
            break;
          }
        }
      }
    }
  }
  return count;
}

int perimeter( const int *ptr_seg, const int &label,
                 const int &len, const int *ptr_nidx ) {
  int count = 0, nidx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] != label ) {
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_nidx[ 6 * i + j ];
        if( nidx == NA_INTEGER ) {
          continue;
        } else {
          if( ptr_seg[ 2 * ( nidx - 1 ) ] == label ) {
            ++ count;
            break;
          }
        }
      }
    }
  }
  return count;
}