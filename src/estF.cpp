#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <set>

#include "helper.h"
#include "initParm.h"
#include "updateBeta.h"
#include "debug.h"
#include "cmpE.h"
#include "copyParm.h"
#include "restoreImg.h"

#include "cnctRegion.h"
#include "excldVoxel.h"
#include "extRegion.h"
#include "pad2zero.h"
#include "regions.h"
#include "region2slice.h"
#include "enclose.h"
#include "inRegion.h"
#include "excldRegion.h"
#include "restoreImg.h"
#include "tissueType.h"
#include "wrapUp.h"
#include "furtherSeg.h"

using std::vector;
using std::set;

// estimate parameters of t1ce and t2 images without tumor
extern "C" SEXP estF( SEXP model, SEXP delta, SEXP gamma, 
                        SEXP alpha, SEXP beta, SEXP lambda2, 
                        SEXP m, SEXP nu2, SEXP maxit );
SEXP estF( SEXP model, SEXP delta, SEXP gamma, 
           SEXP alpha, SEXP beta, SEXP lambda2, 
           SEXP m, SEXP nu2, SEXP maxit ) {
  SEXP info = getListElement( model, "info" );
  SEXP seg = getListElement( model, "seg" );
  SEXP dim_ = getListElement( model, "dim" );
  
  const int *dim__ = INTEGER( dim_ );
  const int *old_seg = INTEGER( seg );
  
  SEXP idx = getListElement( info, "idx" );
  SEXP nidx = getListElement( info, "nidx" );
  SEXP intst = getListElement( info, "intst" );
  SEXP nintst = getListElement( info, "nintst" );
  SEXP aidx = getListElement( info, "aidx" );
  SEXP r_nr = getListElement( info, "nr" );
  SEXP r_nc = getListElement( info, "nc" );
  SEXP r_ns = getListElement( info, "ns" );
  
  const int nr = INTEGER( r_nr )[ 0 ];
  const int nc = INTEGER( r_nc )[ 0 ];
  const int ns = INTEGER( r_ns )[ 0 ];
  
  
  const int *ptr_idx = INTEGER( idx );
  const int *ptr_nidx = INTEGER( nidx );
  const int *ptr_aidx = INTEGER( aidx );
  const double *ptr_intst = REAL( intst );
  const double *ptr_nintst = REAL( nintst );
  
  
  const double *ptr_delta = REAL( delta );
  const double *ptr_gamma = REAL( gamma );
  const double *ptr_alpha = REAL( alpha );
  const double *old_beta = REAL( beta );
  const double *ptr_lambda2 = REAL( lambda2 );
  const double *ptr_m = REAL( m );
  const double *ptr_nu2 = REAL( nu2 );
  const int *ptr_maxit = INTEGER( maxit );
  
  int len = length( idx );
  // number of types of healthy regions
  int n_health = length( alpha );
  
  // a copy of seg and beta
  SEXP res_seg = PROTECT( allocMatrix( INTSXP, 2, len ) );
  SEXP res_beta = PROTECT( allocVector( REALSXP, n_health ) );
  int *ptr_res_seg = INTEGER( res_seg );
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_res_seg[ i ] = old_seg[ i ];
  }
  double *ptr_res_beta = REAL( res_beta );
  for( int i = 0; i < n_health; ++ i ) {
    ptr_res_beta[ i ] = old_beta[ i ];
  }
  // 0, 1, 2 = -1, -2, -3
  // 0, 1 = -1, -2
  vector<double> health_parm( 8 * n_health, 0 );
  
  vector<int> region;
  region.reserve( len );
  vector<double> theta( 6, 0 );
  
  // Initialize parameters for healthy regions
  initParmHealth( n_health, region, theta, health_parm, ptr_res_seg, ptr_m,
                  ptr_nu2, ptr_intst, ptr_lambda2, ptr_nidx,
                  ptr_nintst, ptr_alpha, 
                  ptr_res_beta, len, 20 );
                   
  updateBeta( n_health, ptr_res_beta, ptr_alpha, health_parm );
  
  for( int i = 0; i < *ptr_maxit; ++ i ) {
    for( int j = 1; j <= len; ++ j ) {
      cmpE( n_health, 
            j, health_parm, ptr_res_seg,ptr_nidx, ptr_intst, ptr_nintst, 
            ptr_delta, ptr_gamma, theta ); 
     
    }
    initParmHealth( n_health, region, theta, health_parm, ptr_res_seg, 
                    ptr_m, ptr_nu2, ptr_intst, 
                    ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, 
                    ptr_res_beta, len, 20 );
                     
  }
  // segment zero blocks
  for( int j = 1; j <= len; j ++ ) {
    if( ptr_res_seg[ 2 * ( j - 1 ) ] == 0 ) {
      cmpE( n_health,
            j, health_parm, ptr_res_seg,ptr_nidx, ptr_intst, ptr_nintst, 
            ptr_delta, ptr_gamma, theta );
    }
  }
  
  int nrow = 1 + 2 + 6;
  int ncol = health_parm.size() / 8;
  
  SEXP res_parm = PROTECT( allocMatrix( REALSXP, nrow, ncol ) );
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  double *ptr_res_parm = REAL( res_parm );
  int *ptr_res_image = INTEGER( res_image );
  copyParmHealth( n_health, health_parm, ptr_res_parm, nrow );
  
  // Identify edema and non-enhancing tumor core
  int *ptr_tissue1 = new int[ 2 * len ]();
  int *ptr_tissue2 = new int[ 2 * len ]();
  int *ptr_enclose_ts1 = new int[ 2 * len ]();
  int *ptr_enclose_ts2 = new int[ 2 * len ]();
  int *ptr_seg2_copy = new int[ 2 * len ]();
  int n_ts1 = 0, n_ts2 = 0, n_en_ts1 = 0, n_en_ts2 = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_res_seg[ 2 * i ] == -1 ) {
      ptr_tissue1[ 2 * i ] = 2;
      ++ n_ts1;
    } else if( ptr_res_seg[ 2 * i ] == -2 ){
      ptr_tissue2[ 2 * i ] = 1;
      ++ n_ts2;
    }
  }
  inRegion( ptr_enclose_ts1, len, ptr_tissue2, 1, 
            ptr_tissue1, 2, ptr_seg2_copy, 
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_ts1[ 2 * i ] == 1 ) {
      ++ n_en_ts1;
    }
  }
  inRegion( ptr_enclose_ts2, len, ptr_tissue1, 2,
            ptr_tissue2, 1, ptr_seg2_copy,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_ts2[ 2 * i ] == 1 ) {
      ++ n_en_ts2;
    }
  }
  if( ( n_en_ts1 / ( double )n_ts1 ) > 
      ( n_en_ts2 / ( double )n_ts2 ) ) {
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enclose_ts2[ 2 * i ] == 1 ||
          ptr_tissue1[ 2 * i ] == 2 ) {
        ptr_res_seg[ 2 * i ] = 1;
      } else {
        ptr_res_seg[ 2 * i ] = 2;
      }
      ptr_res_parm[ 0 ] = 1;
      ptr_res_parm[ 9 ] = 2;
    }
  } else {
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enclose_ts1[ 2 * i ] == 1 ||
          ptr_tissue2[ 2 * i ] == 1 ) {
        ptr_res_seg[ 2 * i ] = 1;
      } else {
        ptr_res_seg[ 2 * i ] = 2;
      }
      ptr_res_parm[ 0 ] = 2;
      ptr_res_parm[ 9 ] = 1;
    }
  }           
  
  restoreImg( ptr_idx, ptr_res_seg, ptr_res_image, len, dim__ );
  
  // results to list
  SEXP names = PROTECT( allocVector( STRSXP, 4 ) );
  
  SET_STRING_ELT( names, 0, mkChar( "beta" ) );
  SET_STRING_ELT( names, 1, mkChar( "seg" ) );
  SET_STRING_ELT( names, 2, mkChar( "parm" ) );
  SET_STRING_ELT( names, 3, mkChar( "image" ) );
  
  SEXP res = PROTECT( allocVector( VECSXP, 4 ) );
  SET_VECTOR_ELT( res, 0, res_beta );
  SET_VECTOR_ELT( res, 1, res_seg );
  SET_VECTOR_ELT( res, 2, res_parm );
  SET_VECTOR_ELT( res, 3, res_image );
  setAttrib( res, R_NamesSymbol, names );
  delete [] ptr_tissue1;
  delete [] ptr_tissue2;
  delete [] ptr_enclose_ts1;
  delete [] ptr_enclose_ts2;
  delete [] ptr_seg2_copy;
  
  UNPROTECT( 6 );
  return res;
}