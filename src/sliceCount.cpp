#include "sliceCount.h"
// 0 ~ ( nr - 1 ) : x
// nr ~ ( nr + nc - 1 ): y
// ( nr + nc ) ~ ( nr + nc + ns - 1 ): z
void sliceCount( vector<int> &counts, 
                 int *ptr_seg, const int &len, const int &label,
                 const int *ptr_nidx, const int *ptr_aidx,
                 const int &nr, const int &nc, const int &ns ) {
  int r = 0, c = 0, s = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] == label ) {
      r = ptr_aidx[ 3 * i ];
      c = ptr_aidx[ 3 * i + 1 ];
      s = ptr_aidx[ 3 * i + 2 ];
      
      ++ counts[ r - 1 ];
      ++ counts[ nr + c - 1 ];
      ++ counts[ nr + nc + s - 1 ];
      
    }
  }
}
void sliceCount( vector<int> &counts, const vector<int> &region,
                 const int *ptr_nidx, const int *ptr_aidx,
                 const int &nr, const int &nc, const int &ns ) {
  int r = 0, c = 0, s = 0;
  int idx;
  int size = region.size();
  for( int i = 0; i < size; ++ i ) {
    idx = region[ i ];
    if( idx != 0 ) {
      idx = idx - 1;
      r = ptr_aidx[ 3 * idx ];
      c = ptr_aidx[ 3 * idx + 1 ];
      s = ptr_aidx[ 3 * idx + 2 ];
      
      ++ counts[ r - 1 ];
      ++ counts[ nr + c - 1 ];
      ++ counts[ nr + nc + s - 1 ];
      
    }
  }
}