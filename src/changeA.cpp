#include <R.h>
#include <Rinternals.h>

extern "C" SEXP changeA( SEXP img, SEXP label );

SEXP changeA( SEXP img, SEXP label ) {
  
  double *image = REAL( img );
  int *lbl = INTEGER( label );
  
  SEXP dim = getAttrib( img, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  
  SEXP ans = PROTECT( alloc3DArray( REALSXP, nr, nc, ns ) );
  double *ptr_ans = REAL( ans );
  
  for( int i = 0; i < len; ++ i ) {
    if( lbl[ i ] == 0 ) {
      ptr_ans[ i ] = image[ i ];
    } else {
      ptr_ans[ i ] = R_NaN;
    }
  }
  UNPROTECT( 1 );
  return ans;
}
