#include <R.h>
#include <Rinternals.h>

#include "assignIdxIntst.h"

extern "C" SEXP indexMat( SEXP img, SEXP label );

SEXP indexMat( SEXP img, SEXP label ) {
  int n_nbr = 6;
  const double *image = REAL( img );
  const int *lbl = INTEGER( label );
  
  SEXP dim = getAttrib( img, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  int *ptr_fidx = new int[ len ];
  // Get the number of non-zero voxels and index in valid voxel vector
  int n_valid = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ISNAN( image[ i ] ) || lbl[ i ] == NA_INTEGER ) { 
      // change here for simultaneous flair and t2 model
    // if( ISNAN( flair[ i ] ) || ISNAN( t2[ i ] ) ) { ...
      ptr_fidx[ i ] = NA_INTEGER;
    } else {
      ptr_fidx[ i ] = ++ n_valid; 
    }
  }
  // Get index and intensity of valid voxels
  SEXP valid_idx = PROTECT( allocVector( INTSXP, n_valid ) );
  int *ptr_vidx = INTEGER( valid_idx );
  SEXP valid_intst = PROTECT( allocVector( REALSXP, n_valid ) );
  double *ptr_vintst = REAL( valid_intst );
  n_valid = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_fidx[ i ] != NA_INTEGER ) {
      ptr_vidx[ n_valid ] = i + 1;
      ptr_vintst[ n_valid ++ ] = image[ i ];
    }
  }
  // Get index and intensity of neighboring voxels
  SEXP nbr_idx = PROTECT( allocMatrix( INTSXP, n_nbr, n_valid ) );
  int *ptr_idx = INTEGER( nbr_idx );
  SEXP nbr_intst = PROTECT( allocMatrix( REALSXP, n_nbr, n_valid ) );
  double *ptr_intst = REAL( nbr_intst );
  const int vlen = n_valid;
  n_valid = 0;
  int vec_index = 0;

  for( int k = 0; k < ns; ++ k ) {
    for( int j = 0; j < nc; ++ j ) {
      for( int i = 0; i < nr; ++ i ) {
        vec_index = nc * nr * k + nr * j + i + 1;
        if( vec_index == ptr_vidx[ n_valid ] ) {
          assignIdxIntst( i, j, k, nr, nc, ns, n_nbr, n_valid, ptr_fidx, 
                          ptr_idx, ptr_intst, ptr_vintst );
          ++ n_valid;
          if( n_valid == vlen ) {
            break;
          }
        } 
      }
      if( n_valid == vlen ) {
        break;
      }
    }
    if( n_valid == vlen ) {
      break;
    }
  }

  SEXP names = PROTECT( allocVector( STRSXP, 4 ) );
  
  SET_STRING_ELT( names, 0, mkChar( "idx" ) );
  SET_STRING_ELT( names, 1, mkChar( "nidx" ) );
  SET_STRING_ELT( names, 2, mkChar( "intst" ) );
  SET_STRING_ELT( names, 3, mkChar( "nintst" ) );
  
  SEXP res = PROTECT( allocVector( VECSXP, 4 ) );
  SET_VECTOR_ELT( res, 0, valid_idx );
  SET_VECTOR_ELT( res, 1, nbr_idx );
  SET_VECTOR_ELT( res, 2, valid_intst );
  SET_VECTOR_ELT( res, 3, nbr_intst );
  setAttrib( res, R_NamesSymbol, names );
  delete [] ptr_fidx;
  UNPROTECT( 6 ); 
  return res;
}
