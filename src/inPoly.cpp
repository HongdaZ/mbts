#include <vector>

#include "inPoly.h"

using std::vector;
using std::min;
using std::max;

// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.


// a Point is defined by its coordinates {int x, y;}
//===================================================================

typedef struct {int x, y;} Point;
// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
inline int
isLeft( Point P0, Point P1, Point P2 )
{
  return ( (P1.x - P0.x) * (P2.y - P0.y)
             - (P2.x -  P0.x) * (P1.y - P0.y) );
}
//===================================================================


// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int
cn_PnPoly( Point P, Point* V, int n )
{
  int    cn = 0;    // the  crossing number counter
  
  // loop through all edges of the polygon
  for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
    if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
          || ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
      // compute  the actual edge-ray intersect x-coordinate
      float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
      if (P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
        ++cn;   // a valid crossing of y=P.y right of P.x
    }
  }
  return (cn&1);    // 0 if even (out), and 1 if  odd (in)
  
}
//===================================================================


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int
wn_PnPoly( const Point &P, const vector<Point> &V )
{
  int n = V.size() - 1;
  int    wn = 0;    // the  winding number counter
  
  // loop through all edges of the polygon
  if( V[ n ].x == P.x && V[ n ].y == P.y ) {
    return 3;
  } 
  for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
    if( ( V[ i ].x == P.x && V[ i ].y == P.y ) ||
        ( V[ i + 1 ].x == P.x && V[ i + 1 ].y == P.y )) {
      return 3;
    } else if( ( V[ i ].x - P.x ) * ( V[ i + 1 ].y - P.y  ) ==
      ( V[ i + 1 ].x - P.x ) * ( V[ i ].y - P.y  ) &&
      min( V[ i ].x, V[ i + 1 ].x ) <= P.x && 
      P.x <= max( V[ i ].x, V[ i + 1 ].x ) &&
      min( V[ i ].y, V[ i + 1 ].y ) <= P.y &&
      P.y <= max( V[ i ].y, V[ i + 1 ].y ) ) {
      return 2;
    }else {
      if (V[i].y <= P.y) {          // start y <= P.y
        if (V[i+1].y  > P.y)      // an upward crossing
          if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
            ++wn;            // have  a valid up intersect
      }
      else {                        // start y > P.y (no test needed)
        if (V[i+1].y  <= P.y)     // a downward crossing
          if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
            --wn;            // have  a valid down intersect
      }
    }
  }
  return wn;
}
//===================================================================
vector<int> inPoly( const vector<int> &p, const vector<int> &poly ) {
  int a = 0;
  int len = poly.size() / 2;
  int n_point = p.size() / 2;
  vector<Point> polygon( len );
  vector<int> inside( n_point );
  Point point;
  for( int i = 0; i < len; ++ i ) {
    polygon[ i ].x = poly[ 2 * i ];
    polygon[ i ].y = poly[ 2 * i + 1 ];
  }
  for( int i = 0; i < n_point; ++ i ) {
    point.x = p[ 2 * i ];
    point.y = p[ 2 * i + 1 ];
    a = wn_PnPoly( point, polygon );
    if( a != 2 && a != 3 && a != 0 ) {
      a = 1;
    }
    inside[ i ] = a;
  }
  return inside;
}