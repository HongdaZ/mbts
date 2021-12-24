#ifndef POLYGONAREA_H
#define POLYGONAREA_H

#include <cmath>
#include <vector>

using std::vector;
using std::abs;

// Find the area of polygon
template< typename T >
double polygonArea( const vector<T> &hull ) {
  double area = 0.0;
  const int n = hull.size() / 2;
  int j = n - 1;
  for( int i = 0; i < n; ++ i ) {
    area += ( hull[ 2 * j ] + hull[ 2 * i ] ) * 
      ( hull[ 2 * j + 1 ] - hull[ 2 * i + 1 ] );
    j = i;
  }
  return abs( area / 2.0 );
}

#endif