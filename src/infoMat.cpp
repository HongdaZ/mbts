#include <R.h>
#include <Rinternals.h>

#include "assignIdxIntst.h"

extern "C" SEXP infoMat( SEXP img,
                         SEXP label );
SEXP infoMat( SEXP img,
              SEXP label ) {
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
  SEXP a_idx = PROTECT( allocMatrix( INTSXP, 3, n_valid ) );
  
  int *ptr_idx = INTEGER( nbr_idx );
  int *ptr_aidx = INTEGER( a_idx );
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
          ptr_aidx[ 3 * n_valid ] = i + 1;
          ptr_aidx[ 3 * n_valid + 1 ] = j + 1;
          ptr_aidx[ 3 * n_valid + 2 ] = k + 1;
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
  
  SEXP names = PROTECT( allocVector( STRSXP, 8 ) );
  SEXP r_nr = PROTECT( ScalarInteger( nr ) );
  SEXP r_nc = PROTECT( ScalarInteger( nc ) );
  SEXP r_ns = PROTECT( ScalarInteger( ns ) );
  
  SET_STRING_ELT( names, 0, mkChar( "idx" ) );
  SET_STRING_ELT( names, 1, mkChar( "nidx" ) );
  SET_STRING_ELT( names, 2, mkChar( "intst" ) );
  SET_STRING_ELT( names, 3, mkChar( "nintst" ) );
  SET_STRING_ELT( names, 4, mkChar( "aidx" ) );
  SET_STRING_ELT( names, 5, mkChar( "nr" ) );
  SET_STRING_ELT( names, 6, mkChar( "nc" ) );
  SET_STRING_ELT( names, 7, mkChar( "ns" ) );
  
  SEXP res = PROTECT( allocVector( VECSXP, 8 ) );
  SET_VECTOR_ELT( res, 0, valid_idx );
  SET_VECTOR_ELT( res, 1, nbr_idx );
  SET_VECTOR_ELT( res, 2, valid_intst );
  SET_VECTOR_ELT( res, 3, nbr_intst );
  SET_VECTOR_ELT( res, 4, a_idx );
  SET_VECTOR_ELT( res, 5, r_nr );
  SET_VECTOR_ELT( res, 6, r_nc );
  SET_VECTOR_ELT( res, 7, r_ns );
  setAttrib( res, R_NamesSymbol, names );
  delete [] ptr_fidx;
  UNPROTECT( 10 ); 
  return res;
}