#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "cHull.h"

using std::vector;

// Find the convex hull of a set of points
extern "C" SEXP convexHull( SEXP points );

SEXP convexHull( SEXP points ) {
  const int *ptr_points = INTEGER( points );
  int len = length( points );

  vector<int> a( len );
  for( int i = 0; i < len; ++ i ) {
    a[ i ] = ptr_points[ i ];
  }
  vector<int> res = cHull( a );
  len = res.size() / 2;
  SEXP hull = PROTECT( allocMatrix( INTSXP, 2, len ) );
  int *ptr_hull = INTEGER( hull );
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_hull[ i ] = res[ i ];
  }
  UNPROTECT( 1 );
  return hull;
}