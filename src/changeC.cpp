#include <R.h>
#include <Rinternals.h>

extern "C" SEXP changeC( SEXP old, SEXP image, SEXP q );

SEXP changeC( SEXP old, SEXP image, SEXP q ) {
  int *lbl = INTEGER( old );
  double *img = REAL( image );
  double q_ = REAL( q )[ 0 ];

  SEXP dim = getAttrib( old, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  
  SEXP ans = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  int *ptr_ans = INTEGER( ans );
  
  for( int i = 0; i < len; ++ i ) {
    if( img[ i ] < q_ ) {
      ptr_ans[ i ] = NA_INTEGER;
    } else {
      ptr_ans[ i ] = lbl[ i ];
    }
  }
  UNPROTECT( 1 );
  return ans;
}
