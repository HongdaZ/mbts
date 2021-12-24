#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "inPoly.h"

using std::vector;

// Check if the points are inside the convex hull
// works for closed polygon
extern "C" SEXP inSide( SEXP p, SEXP poly );

SEXP inSide( SEXP p, SEXP poly ) {
  const int *ptr_p = INTEGER( p );
  int len_point = length( p );

  vector<int> a( len_point );
  for( int i = 0; i < len_point; ++ i ) {
    a[ i ] = ptr_p[ i ];
  }
  const int *ptr_poly = INTEGER( poly );
  int len_poly = length( poly );
  
  vector<int> b( len_poly );
  for( int i = 0; i < len_poly; ++ i ) {
    b[ i ] = ptr_poly[ i ];
  }
  
  vector<int> res = inPoly( a, b );
  int len = res.size();
  SEXP inside = PROTECT( allocVector( INTSXP, len ) );
  int *ptr_inside = INTEGER( inside );
  for( int i = 0; i < len; ++ i ) {
    ptr_inside[ i ] = res[ i ];
  }
  UNPROTECT( 1 );
  return inside;
}