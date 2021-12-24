#include <R.h>
#include <Rinternals.h>

extern "C" SEXP call_interface( SEXP x );

SEXP call_interface( SEXP x ) {
  double *ptr_x = REAL( x );
  SEXP ans = PROTECT( allocVector( REALSXP, 1 ) );
  REAL( ans )[ 0 ] = *ptr_x + 1; 
  UNPROTECT( 1 );
  return ans;
}
