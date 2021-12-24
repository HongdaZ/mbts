#include "cnctRegion.h"
#include "findRegion.h"

// Find the connected region or return false
// start is an R style index (start = i + 1) 
bool cnctRegion( const int &start, const int *ptr_nidx,
                 int *ptr_seg1, int *ptr_seg2, 
                 const int &label, vector<int> &region ) {
  bool res;
  if( ptr_seg1[ 2 * ( start - 1 ) ] > 0 && 
      ptr_seg1[ 2 * start - 1 ] == 0 && 
      ptr_seg2[ 2 * ( start - 1 ) ] == label ) {
    res = true;
    int n_region = 1;
    list<int> tumor_nbr;
    bool early_return;
    findRegion( region, n_region, ptr_seg2, ptr_nidx, true, 
                tumor_nbr, early_return, start );
    for( vector<int>::iterator it = region.begin(); it != region.end();
    ++ it ) {
      ptr_seg1[ 2 * *it - 1 ] = 1;
    }
  } else {
    ptr_seg1[ 2 * start - 1 ] = 1;
    res = false;
  }
  return res;
}
// 2D version of function above
bool cnctRegion( const int &start, const int *ptr_nidx, 
                 const int *ptr_aidx, const int &plane,
                 int *ptr_seg1, int *ptr_seg2, 
                 const int &label, vector<int> &region ) {
  bool res;
  if( ptr_seg1[ 2 * ( start - 1 ) ] > 0 && 
      ptr_seg1[ 2 * start - 1 ] == 0 && 
      ptr_seg2[ 2 * ( start - 1 ) ] == label ) {
    res = true;
    int n_region = 1;
    list<int> tumor_nbr;
    bool early_return;
    findRegion( region, n_region, ptr_seg2, ptr_nidx, ptr_aidx, plane,
                true, tumor_nbr, early_return, start );
    for( vector<int>::iterator it = region.begin(); it != region.end();
    ++ it ) {
      ptr_seg1[ 2 * *it - 1 ] = 1;
    }
  } else {
    ptr_seg1[ 2 * start - 1 ] = 1;
    res = false;
  }
  return res;
}