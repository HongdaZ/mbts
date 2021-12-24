#include <cmath>

#include "roundness.h"
#include "radius.h"

using std::sqrt;

// Measure the roundness of the voxels in region
double roundness( const vector<int> &region, const int *ptr_aidx ) {
  vector<double> shift{ 0.5, 0.5, 0.5,
                        0.5, 0.5, -0.5,
                        0.5, -0.5, 0.5,
                        0.5, -0.5, -0.5,
                        -0.5, 0.5, 0.5,
                        -0.5, 0.5, -0.5,
                        -0.5, -0.5, 0.5,
                        -0.5, -0.5, -0.5 };
  double r = radius( region.size() );
  double max_dist = 0, dist = 0;
  const int len = region.size();
  double cr = 0, cc = 0, cs = 0, r1, c1, s1, n = 0;
  double r2, c2, s2;
  int index1;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
      c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
      s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
      for( int j = 0; j < 8; ++ j ) {
        cr += r1 + shift[ 3 * j ];
        cc += c1 + shift[ 3 * j + 1 ];
        cs += s1 + shift[ 3 * j + 2 ];
        ++ n;
      }
    }
  }
  cr /= n;
  cc /= n;
  cs /= n;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
      c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
      s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
      for( int j = 0; j < 8; ++ j ) {
        r2 = r1 + shift[ 3 * j ];
        c2 = c1 + shift[ 3 * j + 1 ];
        s2 = s1 + shift[ 3 * j + 2 ];
        dist =  sqrt( pow( r2 - cr, 2 ) + pow( c2 - cc, 2 ) +
          pow( s2 - cs, 2 ) );
        if( max_dist < dist ) {
          max_dist = dist;
        }
      }
    }
  }
  
  max_dist /= r;
  return max_dist * max_dist;
}
// 2D version of the function above
double roundness( const vector<int> &region, const int &plane,
               const int *ptr_aidx ) {
  vector<double> shift{ 0.5, 0.5,
                        0.5, -0.5,
                        -0.5, 0.5,
                        -0.5, -0.5 };
  vector<int> plane_idx{ 1, 2, 0, 2, 0, 1 };
  vector<int> curr_plane{ plane_idx[ 2 * plane ], 
                          plane_idx[ 2 * plane + 1 ] };
  double r = radius2D( region.size() );
  double max_dist = 0, dist = 0;
  const int len = region.size();
  double n = 0;
  vector<double> center( 2, 0 );
  vector<double> xyz( 2, 0 );
  int index1;
  double tmp1, tmp2;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      for( int j = 0; j < 2; ++ j ) {
        xyz[ j ] = ptr_aidx[ 3 * ( index1 - 1 ) + curr_plane[ j ] ];
        for( int k = 0; k < 4; ++ k ) {
          center[ j ] += xyz[ j ] + shift[ 2 * k + j ];
          ++ n;
        }
      }
    }
  }
  n /= 2;
  center[ 0 ] /= n;
  center[ 1 ] /= n;
  
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      xyz[ 0 ] = ptr_aidx[ 3 * ( index1 - 1 ) + curr_plane[ 0 ] ];
      xyz[ 1 ] = ptr_aidx[ 3 * ( index1 - 1 ) + curr_plane[ 1 ] ];
      for( int j = 0; j < 4; ++ j ) {
        tmp1 = xyz[ 0 ] + shift[ 2 * j ];
        tmp2 = xyz[ 1 ] + shift[ 2 * j + 1 ];
        dist =  sqrt( pow( tmp1 - center[ 0 ], 2 ) + 
          pow( tmp2 - center[ 1 ], 2 ) );
        if( max_dist < dist ) {
          max_dist = dist;
        }
      }
    }
  }
  max_dist /= r;
  return max_dist * max_dist;
}