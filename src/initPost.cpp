#include <R.h>
#include <Rinternals.h>

#include "assignIdx.h"

extern "C" SEXP initPost( SEXP t1ce_image, SEXP flair_image, 
                         SEXP t2_image, SEXP t1ce_intst );
SEXP initPost( SEXP t1ce_image, SEXP flair_image, 
               SEXP t2_image, SEXP t1ce_intst ) {
  int n_nbr = 6;
  const int *ptr_t1ce = INTEGER( t1ce_image );
  const int *ptr_flair = INTEGER( flair_image );
  const int *ptr_t2 = INTEGER( t2_image );
  const double *ptr_t1ce_intst = REAL( t1ce_intst );
  
  SEXP dim = getAttrib( t1ce_image, R_DimSymbol );
  int nr = INTEGER( dim )[ 0 ];
  int nc = INTEGER( dim )[ 1 ];
  int ns = INTEGER( dim )[ 2 ];
  int len = nr * nc * ns;
  
  // Get the number of non-zero voxels and indices
  int *fidx = new int[ len ];
  int n_valid = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_t1ce[ i ] == 0 && ptr_flair[ i ] == 0 &&
        ptr_t2[ i ] == 0 ) {
      fidx[ i ] = NA_INTEGER;
    } else {
      fidx[ i ] = ++ n_valid;
    }
  }
  
  // Get indices and labels of valid voxels
  SEXP vidx = PROTECT( allocVector( INTSXP, n_valid ) );
  SEXP t1ce_seg = PROTECT( allocMatrix( INTSXP, 2, n_valid ) );
  SEXP flair_seg = PROTECT( allocMatrix( INTSXP, 2, n_valid ) );
  SEXP t2_seg = PROTECT( allocMatrix( INTSXP, 2, n_valid ) );
  SEXP t1ce_grey = PROTECT( allocVector( REALSXP, n_valid ) );
  int *ptr_vidx = INTEGER( vidx );
  int *ptr_t1ce_seg = INTEGER( t1ce_seg );
  int *ptr_flair_seg = INTEGER( flair_seg );
  int *ptr_t2_seg = INTEGER( t2_seg );
  double *ptr_t1ce_grey = REAL( t1ce_grey );
  n_valid = 0;
  for( int i = 0; i < len; ++ i ) {
    if( fidx[ i ] != NA_INTEGER ) {
      ptr_vidx[ n_valid ] = i + 1;
      ptr_t1ce_seg[ n_valid * 2 ] = ptr_t1ce[ i ];
      ptr_flair_seg[ n_valid * 2] = ptr_flair[ i ];
      ptr_t2_seg[ n_valid * 2 ] = ptr_t2[ i ];
      ptr_t1ce_grey[ n_valid ] = ptr_t1ce_intst[ i ];
      ++ n_valid;
    }
  }
  
  // Get array index and indices of neighboring voxels
  SEXP nbr_idx = PROTECT( allocMatrix( INTSXP, n_nbr, n_valid ) );
  SEXP a_idx = PROTECT( allocMatrix( INTSXP, 3, n_valid ) );
  
  int *ptr_nidx = INTEGER( nbr_idx );
  int *ptr_aidx = INTEGER( a_idx );
  const int vlen = n_valid;
  n_valid = 0;
  int vec_index = 0;
  
  for( int k = 0; k < ns; ++ k ) {
    for( int j = 0; j < nc; ++ j ) {
      for( int i = 0; i < nr; ++ i ) {
        vec_index = nc * nr * k + nr * j + i + 1;
        if( vec_index == ptr_vidx[ n_valid ] ) {
          // Assign array index
          ptr_aidx[ 3 * n_valid ] = i + 1;
          ptr_aidx[ 3 * n_valid + 1 ] = j + 1;
          ptr_aidx[ 3 * n_valid + 2 ] = k + 1;
          assignIdx( i, j, k, nr, nc, ns, n_nbr, n_valid, fidx, 
                     ptr_nidx );
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
  
  SEXP names = PROTECT( allocVector( STRSXP, 10 ) );
  SEXP r_nr = PROTECT( ScalarInteger( nr ) );
  SEXP r_nc = PROTECT( ScalarInteger( nc ) );
  SEXP r_ns = PROTECT( ScalarInteger( ns ) );
  
  SET_STRING_ELT( names, 0, mkChar( "t1ce_seg" ) );
  SET_STRING_ELT( names, 1, mkChar( "flair_seg" ) );
  SET_STRING_ELT( names, 2, mkChar( "t2_seg" ) );
  SET_STRING_ELT( names, 3, mkChar( "idx" ) );
  SET_STRING_ELT( names, 4, mkChar( "nidx" ) );
  SET_STRING_ELT( names, 5, mkChar( "aidx" ) );
  SET_STRING_ELT( names, 6, mkChar( "nr" ) );
  SET_STRING_ELT( names, 7, mkChar( "nc" ) );
  SET_STRING_ELT( names, 8, mkChar( "ns" ) );
  SET_STRING_ELT( names, 9, mkChar( "t1ce_intst" ) );
  
  SEXP res = PROTECT( allocVector( VECSXP, 10 ) );
  SET_VECTOR_ELT( res, 0, t1ce_seg );
  SET_VECTOR_ELT( res, 1, flair_seg );
  SET_VECTOR_ELT( res, 2, t2_seg );
  SET_VECTOR_ELT( res, 3, vidx );
  SET_VECTOR_ELT( res, 4, nbr_idx );
  SET_VECTOR_ELT( res, 5, a_idx );
  SET_VECTOR_ELT( res, 6, r_nr );
  SET_VECTOR_ELT( res, 7, r_nc );
  SET_VECTOR_ELT( res, 8, r_ns );
  SET_VECTOR_ELT( res, 9, t1ce_grey );
  setAttrib( res, R_NamesSymbol, names );
  delete [] fidx;
  UNPROTECT( 12 ); 
  return res;
}