#include "pTNbr.h"

#include "Rinternals.h"

// Find the proportion of new voxels having tumor neighbor
double pTNbr( const vector<int> &region, const int *ptr_seg1, 
           const int &label1, const int *ptr_nidx ) {
  int total = 0, n = 0, idx = 0, n_idx = 0;
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ++ total;
      for( int i = 0; i < 6; ++ i ) {
        n_idx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( n_idx !=  NA_INTEGER ) {
          if( ptr_seg1[ 2 * ( n_idx - 1 ) ] == label1 ) {
            ++ n;
            break;
          }
        }
      } 
    }
  }
  return (double)n / total;
}
// 2D version of the function above
double pTNbr2D( const vector<int> &region, const int *ptr_seg1, 
                const int &label1, const int *ptr_nidx, 
                const int *ptr_aidx, const int &plane ) {
  int total = 0, n = 0, idx = 0, n_idx = 0, slice_region, slice_nbr;
  vector<int>::const_iterator it = region.begin();
  idx = *it;
  slice_region = ptr_aidx[ 3 * ( idx - 1 ) + plane ];
  for( ; it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ++ total;
      for( int i = 0; i < 6; ++ i ) {
        n_idx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( n_idx !=  NA_INTEGER ) {
          if( ptr_seg1[ 2 * ( n_idx - 1 ) ] == label1 ) {
            slice_nbr = ptr_aidx[ 3 * ( n_idx - 1 ) + plane ]; 
            if( slice_nbr == slice_region ) {
              ++ n;
              break;
            }
          }
        }
      } 
    }
  }
  return (double)n / total;
}