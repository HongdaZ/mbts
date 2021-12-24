#include "regions.h"
#include "cnctRegion.h"
#include "pad2zero.h"

// Find all the connected regions in ptr_seg with label == label
list<vector<int>> regions( int *ptr_seg, const int &len, 
                           vector<int> &region,
                           const int &label,
                           const int *ptr_nidx,
                           const int *ptr_aidx ) {
  list<vector<int>> res;
  int n;
  int idx;
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_seg, ptr_seg, 
                    label, region ) ) {
      n = region.size();
      vector<int> r( n * 4 );
      for( int j = 0; j < n; ++ j ) {
        idx = region[ j ];
        r[ 4 * j ] = idx;
        for( int k = 0; k < 3; ++ k ) {
          r[ 4 * j + k + 1 ] = ptr_aidx[ 3 * ( idx - 1 ) + k ];
        }
      }
      res.push_back( r );
    }
  }
  pad2zero( ptr_seg, len );
  return res;
}
// 2D version of function above
vector<list<vector<int>>> regions2D( int *ptr_seg, const int &len, 
                                     vector<int> &region,
                                     const int &label,
                                     const int *ptr_nidx,
                                     const int *ptr_aidx ) {
  vector<list<vector<int>>> res( 3 );
  int n;
  int idx;
  for( int k = 0; k < 3; ++ k ) {
    list<vector<int>> plane;
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, k, ptr_seg, ptr_seg, 
                      label, region ) ) { 
        n = region.size();
        vector<int> r( n * 4 );
        for( int j = 0; j < n; ++ j ) {
          idx = region[ j ];
          r[ 4 * j ] = idx;
          for( int k = 0; k < 3; ++ k ) {
            r[ 4 * j + k + 1 ] = ptr_aidx[ 3 * ( idx - 1 ) + k ];
          }
        }
        plane.push_back( r );
      }
    }
    res[ k ] = plane;
    pad2zero( ptr_seg, len );
  }
  return res;
}
// All voxels with label == label is a whole region
list<vector<int>> regions( const int *ptr_seg, const int &len, 
                           const int &label,
                           const int *ptr_aidx ) {
  list<vector<int>> res;
  list<int> tissue;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] == label ) {
      tissue.push_back( i + 1 );
    }
  }
  const int n = tissue.size();
  vector<int> r( n * 4 );
  int i = 0;
  int idx = 0;
  for( list<int>::const_iterator it_t = tissue.begin();
       it_t != tissue.end(); ++ it_t, ++ i  ) {
    idx = *it_t;
    r[ 4 * i ] = idx;
    for( int j = 0; j < 3; ++ j ) {
      r[ 4 * i + 1 + j ] = ptr_aidx[ 3 * ( idx - 1 ) + j ];
    }
  }
  res.push_back( r );
  return res;
}