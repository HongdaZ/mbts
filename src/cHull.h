#ifndef CHULL_H
#define CHULL_H

#include <vector>
#include <cassert>
#include "quickhull.hpp"
#include "point.hpp"

using std::vector;

// Find convex hull of points a
// Find convex hull
template< typename T >
vector<T> cHull( const vector<T> &a ) {
  using P = point<T>;
  using S = std::vector<P>;
  using I = typename S::iterator;
  
  int len = a.size() / 2;
  S bag( len );
  
  for( int i = 0; i < len; ++ i ) {
    bag[ i ].x = a[ 2 * i ];
    bag[ i ].y = a[ 2 * i + 1 ];
  }
  I rest = quickhull::solve( bag.begin(), bag.end() );
  auto h = rest - bag.begin();
  assert( h == 2 );
  
  vector<T> res( 2 * h );
  for( int i = 0; i < h; ++ i ) {
    res[ 2 * i ] = bag[ i ].x;
    res[ 2 * i + 1 ] = bag[ i ].y;
  }
  return res;
}

#endif