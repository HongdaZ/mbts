#include <R.h>
#include <Rinternals.h>

#include <cmath>   

#include "updateParm.h"
#include "labelRegion.h" 
#include "initMV.h"
#include "zeroVector.h"

using std::abs;

// previous updateParm is slow because of using std::find with std::list
// and extra times of copying matrices

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const vector<int> &region,
                 const double &m,
                 const double &nu2,
                 const double *ptr_intst,
                 const int &curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double &alphal,
                 const double &betal,
                 const int &maxit ) {

  mu = -1;
  zeroVector( theta );
  int i = 0;
  double tol = 1;
  double tmp = 0;
  // Initialize matrices and vectors
  int len_region = region.size();
  int nrow = len_region;
  int ncol = 6;
  double *yln_ = new double[ nrow * ncol ];
  double *yln_i = new double[ nrow * ncol ];
  double *yl_ = new double[ nrow ];
  double sum_y = 0;
  
  initMV( region, yln_, yln_i,yl_, sum_y, ptr_intst, ptr_nidx, ptr_nintst, 
          ptr_seg, curr_label );
  
  double *yln = new double[ nrow * ncol ];
  double *ylna = new double[ nrow * ( ncol / 2 ) ];
  double *yl = new double[ nrow ];
  
  double sum_theta;
  while( i < 1 || ( i < maxit && tol > .0001 )  ) {
    sum_theta = 0;
    for( int j = 0; j < 6; ++ j ) {
      sum_theta += theta[ j ];
    }
    tmp = updateMu( nrow, sigma2, m, nu2, sum_theta, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( nrow, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal, yl, yln, ylna,
              yln_, yln_i, yl_ );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  delete [] yln_;
  delete [] yln_i;
  delete [] yl_;
  delete [] yln;
  delete [] ylna;
  delete [] yl;
  
  return;
}

// update parameters for tumor cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const vector<int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 const int &curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 const int &maxit ) {
  mu = -1;
  zeroVector( theta );
  int i = 0;
  double tol = 1;
  double tmp = 0;
  int len_region = region.size();
  // the voxel in region is labelled as 3 in ptr_seg[ 2, ]
  labelRegion( region, ptr_seg );
  
  // Initialize matrices and vectors
  int nrow = len_region;
  int ncol = 6;
  double *yln_ = new double[ nrow * ncol ];
  double *yln_i = new double[ nrow * ncol ];
  double *yl_ = new double[ nrow ];
  
  double sum_y = 0;
  
  initMV( region, yln_, yln_i,yl_, sum_y, ptr_intst, ptr_nidx, ptr_nintst, 
          ptr_seg, curr_label );
  
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  double *ylna = new double[ nrow * ( ncol / 2 ) ];
  double sum_theta;
  
  while(  i < 1 || ( i < maxit && tol > .0001 ) ) {
    sum_theta = 0;
    for( int j = 0; j < 6; ++ j ) {
      sum_theta += theta[ j ];
    }
    tmp = updateMu( nrow, sigma2, m, mk_1, a, b, sum_theta, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( nrow, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal, yl, yln, ylna,
              yln_, yln_i, yl_ );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  // change ptr_seg[ 2, region ] back
  recoverLabel( region, ptr_seg );
  delete [] yln_;
  delete [] yln_i;
  delete [] yl_;
  delete [] yln;
  delete [] ylna;
  delete [] yl;
  return;
}

// update parameters for outliers or single point tumor regions
void updateParm( double &mu, double &sigma2, 
                 const int idx,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 int *ptr_seg,
                 const double alphal,
                 const double betal,
                 const int &maxit ) {
  mu = -1;
  int i = 0;
  double tol = 1;
  double tmp = 0;
  double sum_y = ptr_intst[ idx - 1 ];
  int n = 1;
  while(  i < 1 || ( i < maxit && tol > .0001 ) ) {
    tmp = updateMu( n, sigma2, m, mk_1, a, b, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateSigma( idx, mu, ptr_intst, alphal, betal, sigma2 );
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  return;
}