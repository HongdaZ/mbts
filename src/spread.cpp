#include "Rinternals.h"

#include <cmath>

#include "spread.h"
#include "radius.h"
#include "perimeter.h"
#include "zeroVector.h"

using std::sqrt;

// Measure the spread of the voxels in region
double spread( const vector<int> &region, int *ptr_seg_copy,
               const int &len, const int *ptr_nidx ) {
  const double PI = 3.14159265358979323846;
  double r = radius( region.size() );
  double spread_idx = 0;
  int idx = 0;
  zeroVector( ptr_seg_copy, len );
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ptr_seg_copy[ 2 * ( idx - 1 ) ] = 1;
    }
  }
  int p = perimeter( ptr_seg_copy, 1, len, ptr_nidx );
  spread_idx = sqrt( ( double ) p / 4 / PI ) / r;
  
  return spread_idx * spread_idx;
}
// 2D version of the function above
double spread( const vector<int> &region, int *ptr_seg_copy,
               const int &len, const int *ptr_nidx,
               const int &plane ) {
  const double PI = 3.14159265358979323846;
  double r = radius2D( region.size() );
  double spread_idx = 0;
  int idx = 0;
  zeroVector( ptr_seg_copy, len );
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ptr_seg_copy[ 2 * ( idx - 1 ) ] = 1;
    }
  }
  int p = perimeter( ptr_seg_copy, 1, len, ptr_nidx, plane );
  spread_idx = ( double ) p / 2 / PI / r;
  return spread_idx * spread_idx;
}