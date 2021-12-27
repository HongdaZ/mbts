#include <R.h>
#include <Rinternals.h>
#include "R_ext/Rdynload.h"
#include "header.h"

extern "C" SEXP call_interface( SEXP x );

attribute_hidden SEXP call_interface( SEXP x ) {
  double *ptr_x = REAL( x );
  SEXP ans = PROTECT( allocVector( REALSXP, 1 ) );
  REAL( ans )[ 0 ] = *ptr_x + 1; 
  UNPROTECT( 1 );
  return ans;
}
