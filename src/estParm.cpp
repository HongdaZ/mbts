#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <chrono> 

#include "helper.h"
#include "initRegion.h"
#include "scTrn.h"
#include "initParm.h"
#include "cmpET.h"
#include "cmpE.h"
#include "updateBeta.h"
#include "skip.h"
#include "debug.h"
#include "copyParm.h"
#include "restoreImg.h"

using std::vector;
using namespace std::chrono; 

// predict the labels for T1ce or T2 images
extern "C" SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
                      SEXP alpha, SEXP beta, SEXP lambda2, 
                      SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit );

SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
            SEXP alpha, SEXP beta, SEXP lambda2, 
            SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit ) {
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
  const double *ptr_a = REAL( a );
  const double *ptr_b = REAL( b );
  const double *ptr_m = REAL( m );
  const double *ptr_nu2 = REAL( nu2 );
  const int *ptr_maxit = INTEGER( maxit );
  
  int len = length( idx );
  // number of types of healthy regions
  int n_health = length( alpha );
  // a copy of seg and beta
  SEXP res_seg = PROTECT( allocMatrix( INTSXP, 2, len ) );
  SEXP res_beta = PROTECT( allocVector( REALSXP, 4 ) );
  int *ptr_res_seg = INTEGER( res_seg );
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_res_seg[ i ] = old_seg[ i ];
  }
  double *ptr_res_beta = REAL( res_beta );
  for( int i = 0; i < 4; ++ i ) {
    ptr_res_beta[ i ] = old_beta[ i ];
  }
  vector<int> tumor_labels( len, 0 );
  vector<int> outl_labels( len, 0 );
  int n_tumor = 0;
  int n_outl = 0;
  // 0, 1, 2 = -1, -2, -3
  vector<double> health_parm( 8 * 3, 0 );
  // 0, 1, 2, ... = -4, -5, -6, ...
  vector<double> tumor_parm( 8 * len, 0 );
  // 0, 1, 2, ... = 1, 2, 3, ...
  vector<double> outl_parm( 2 * len, 0 );
  list<list<int>> tumor_regions;
  
  // whole and subregions in sc functions
  vector<int> regions_whole;
  regions_whole.reserve( len );
  vector<int> regions_sub;
  regions_sub.reserve( len );
  
  // frontier in search
  vector<int> front;
  front.reserve( len );
  // store the result of findRegion
  vector<int> region;
  region.reserve( len );
  
  int n_row = 1 + 2 + 6;
  int n_col = 0;
  
  vector<double> outlier_parm( 3, 0 );
  vector<double> theta( 6, 0 );
  vector<double> tmp_parm( 9, 0 );
  vector<double> out_theta( 6, 0 );
  vector<double> new_out_parm( 2, 0 );
  vector<double> whole_parm( 8, 0 );
  vector<double> label_whole_parm( 9, 0 );
  vector<double> region_parm;
  region_parm.reserve( n_row * 7 );
  
  initRegion( region, tumor_regions,
              ptr_res_seg, ptr_nidx, len,
              tumor_labels, n_tumor );
  initParm( n_health,
            region, theta, true, health_parm, tumor_parm, ptr_res_seg,
            ptr_m, ptr_nu2, ptr_intst, ptr_lambda2, ptr_nidx,
            ptr_nintst, ptr_alpha, ptr_res_beta, tumor_regions,
            ptr_a, ptr_b, len, 20 );
  updateBeta( n_health,
              ptr_res_beta, ptr_alpha, health_parm, tumor_regions,
              tumor_parm );
  
  int old_label = 0;
  int new_label = 0;
  int curr_idx = 0;
  double curr_intst = 0;
  int curr_label = 0;
  int n_region;
  
  auto start = high_resolution_clock::now(); 
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 
  int spend = 0; 
  int max_idx;
  list<int> tumor_nbr;
  list<int> tumor_label;
  int sc;
  bool skip_;
  // whether getParm or updateParm
  list<int> update_parm;
  for( int i = 0; i < *ptr_maxit; ++ i ) {
    for( int j = 1; j <= len; ++ j ) {
      curr_idx = j;
      skip_ = skip( curr_idx, ptr_res_seg, ptr_nidx );
      if( ! skip_ ) {
        curr_label = ptr_res_seg[ 2 * ( curr_idx - 1 ) ];
        if( curr_label < 1 && curr_label > -4 ) {
          cmpE( n_health,
                curr_idx, health_parm, ptr_res_seg, ptr_nidx, ptr_intst,
                ptr_nintst, ptr_delta, ptr_gamma, theta );
        } else {
          sc = scTrn( n_region, update_parm, tumor_regions, region, 
                       tumor_labels, ptr_res_seg, ptr_nidx, 
                       len, curr_idx, regions_whole, regions_sub,
                       tumor_nbr, tumor_label );
          cmpET( n_health, 
                 update_parm,
                 region, curr_idx, sc, regions_whole, regions_sub,
                 tumor_labels, outl_labels, health_parm,
                 tumor_parm, outl_parm, 
                 ptr_res_seg, ptr_nidx, ptr_intst, ptr_nintst, 
                 ptr_delta, ptr_gamma, ptr_alpha, ptr_res_beta, 
                 ptr_lambda2, ptr_a, ptr_b, ptr_m,
                 ptr_nu2, outlier_parm, theta, tmp_parm, out_theta,
                 new_out_parm, whole_parm, label_whole_parm,
                 region_parm, n_tumor, n_outl, 
                 n_region, tumor_regions, n_row, tumor_label );
        }
      }
    }
    Rprintf( "update parm for healthy and tumorous\n" );
    //update parm for healthy and tumorous regions
    initParm( n_health,
              region, theta, false, health_parm, tumor_parm, ptr_res_seg,
              ptr_m, ptr_nu2, ptr_intst, ptr_lambda2, ptr_nidx,
              ptr_nintst, ptr_alpha, ptr_res_beta, tumor_regions,
              ptr_a, ptr_b, len, 20 );
  }
  // segment zero blocks
  for( int j = 1; j <= len; j ++ ) {
    if( ptr_res_seg[ 2 * ( j - 1 ) ] == 0 ) {
      curr_idx = j;
      sc = scTrn( n_region, update_parm, tumor_regions, region,
                   tumor_labels, ptr_res_seg, ptr_nidx,
                   len, curr_idx, regions_whole, regions_sub,
                   tumor_nbr, tumor_label );
      cmpET( n_health,
             update_parm, 
             region, curr_idx, sc, regions_whole, regions_sub,
             tumor_labels, outl_labels, health_parm,
             tumor_parm, outl_parm,
             ptr_res_seg, ptr_nidx, ptr_intst, ptr_nintst,
             ptr_delta, ptr_gamma, ptr_alpha, ptr_res_beta,
             ptr_lambda2, ptr_a, ptr_b, ptr_m,
             ptr_nu2, outlier_parm, theta, tmp_parm, out_theta,
             new_out_parm, whole_parm, label_whole_parm,
             region_parm, n_tumor, n_outl,
             n_region, tumor_regions, n_row, tumor_label );
    }
  }
  n_col = n_tumor + n_health + n_outl;
  
  SEXP res_parm = PROTECT( allocMatrix( REALSXP, n_row, n_col ) );
  SEXP res_image = PROTECT( alloc3DArray( INTSXP,
                                          dim__[ 0 ], 
                                          dim__[ 1 ], 
                                          dim__[ 2 ] ) );
  double *ptr_res_parm = REAL( res_parm );
  int *ptr_res_image = INTEGER( res_image );
  copyParm( n_health,
            health_parm, tumor_parm, outl_parm, ptr_res_parm, n_row,
            tumor_regions, outl_labels, len );
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