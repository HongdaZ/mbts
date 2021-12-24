#include <R.h>
#include <Rinternals.h>

extern "C" SEXP changeB( SEXP label, SEXP image1, SEXP q1, SEXP image2, 
                        SEXP q2, SEXP k  );

SEXP changeB( SEXP label, SEXP image1, SEXP q1, SEXP image2, SEXP q2, 
              SEXP k  ) {
  int *lbl = INTEGER( label );
  double *img1 = REAL( image1 );
  double *img2 = REAL( image2 );
  double q_1 = REAL( q1 )[ 0 ];
  double q_2 = REAL( q2 )[ 0 ];
  int l = INTEGER( k )[ 0 ];
  
  SEXP dim = getAttrib( label, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  
  SEXP ans = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  int *ptr_ans = INTEGER( ans );

  for( int i = 0; i < len; ++ i ) {
    if( img1[ i ] < q_1 && img2[ i ] > q_2 ) {
      ptr_ans[ i ] = l;
    } else {
      ptr_ans[ i ] = lbl[ i ];
    }
  }

  UNPROTECT( 1 );
  return ans;
}
