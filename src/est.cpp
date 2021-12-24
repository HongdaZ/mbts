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

using std::vector;
using std::set;

// estimate parameters of t1ce and t2 images without tumor
extern "C" SEXP est( SEXP model, SEXP delta, SEXP gamma, 
                        SEXP alpha, SEXP beta, SEXP lambda2, 
                        SEXP m, SEXP nu2, SEXP maxit );
SEXP est( SEXP model, SEXP delta, SEXP gamma, 
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
  
  const int *ptr_idx = INTEGER( idx );
  const int *ptr_nidx = INTEGER( nidx );
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
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, 
                                          dim__[ 0 ], 
                                          dim__[ 1 ], 
                                          dim__[ 2 ] ) );
  double *ptr_res_parm = REAL( res_parm );
  int *ptr_res_image = INTEGER( res_image );
  copyParmHealth( n_health, health_parm, ptr_res_parm, nrow );
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
  
  UNPROTECT( 6 );
  return res;
}