#include <cmath>
#include <limits>

#include "center.h"

using std::pow;
using std::sqrt;

// Return the center voxel of outer
int center( const list<int> &outer, const int *ptr_aidx ) {
  double dist = 0;
  double min_dist = std::numeric_limits<double>::infinity();
  double cr = 0, cc = 0, cs = 0, r1, c1, s1;
  int n = outer.size();
  int index1;
  int center;
  list<double> outer_dist;
  for( list<int>::const_iterator it = outer.begin();
       it != outer.end(); ++ it ) {
    index1 = *it;
    r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
    c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
    s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
    cr += r1;
    cc += c1;
    cs += s1;
  }
  cr /= n;
  cc /= n;
  cs /= n;
  for( list<int>::const_iterator it = outer.begin();
       it != outer.end(); ++ it ) {
    index1 = *it;
    r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
    c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
    s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
    dist =  sqrt( pow( r1 - cr + 0.71, 2 ) + pow( c1 - cc + 0.71, 2 ) +
      pow( s1 - cs + 0.71, 2 ) );
    outer_dist.push_back( dist );
    if( dist < min_dist ) {
      min_dist = dist;
    }
  }
  list<int>::const_iterator it_o = outer.begin();
  list<double>::const_iterator it_d = outer_dist.begin();
  for( ; it_d != outer_dist.end(); ++ it_d, ++ it_o ) {
    dist = *it_d;
    if( dist == min_dist ) {
      center = *it_o;
      break;
    }
  }
  return center;
}