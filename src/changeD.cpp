#include <R.h>
#include <Rinternals.h>

extern "C" SEXP changeD( SEXP old, SEXP a, SEXP b );

SEXP changeD( SEXP old, SEXP a, SEXP b ) {
  int *lbl = INTEGER( old );
  int a_ = INTEGER( a )[ 0 ];
  int b_ = INTEGER( b )[ 0 ];
  
  SEXP dim = getAttrib( old, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  
  SEXP ans = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  int *ptr_ans = INTEGER( ans );
  
  for( int i = 0; i < len; ++ i ) {
    if( lbl[ i ] == a_ || lbl[ i ] == b_ ) {
      ptr_ans[ i ] = NA_INTEGER;
    } else {
      ptr_ans[ i ] = lbl[ i ];
    }
  }
  UNPROTECT( 1 );
  return ans;
}
