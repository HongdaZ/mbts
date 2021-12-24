#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <vector>
#include <algorithm>

#include "energyY.h"
#include "labelRegion.h"
#include "initMV.h"

using std::vector;
using std::find;

# define Pi 3.14159265358979323846

// region starts from 1
// calculate energy for tumor regions
double energyY( const vector<int> &region,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                const double &alphak,
                const double &betak,
                const double &a, const double &b ){
  // the voxel in region is labelled as 3 in ptr_seg[ 2, ]
  labelRegion( region, ptr_seg );
  
  int nrow = region.size();
  int ncol = 6;
  
  double *yln_ = new double[ nrow * ncol ];
  double *yln_i = new double[ nrow * ncol ];
  double *yl_ = new double[ nrow ];
  double sum_y = 0;
  // Initialize yln_, yln_i, and yl;
  initMV( region, yln_, yln_i,yl_, sum_y, ptr_intst, ptr_nidx, ptr_nintst,
          ptr_seg, - 4 );
  
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  // Initialize yln and yl;
  initMV( yl, yl_, yln, yln_, yln_i, nrow, mu ); 
  
  // yl - yln * thetal
  double *vtheta = new double[ ncol ];
  for( int i = 0; i < 6; ++ i ) {
    vtheta[ i ] = theta[ i ];
  }
  int incx = 1;
  double alpha = - 1;
  double beta = 1;
  F77_CALL( dgemv )( "n", &nrow, &ncol, &alpha, yln, &nrow, vtheta, &incx,
            &beta, yl, &incx );
  
  double energy = 0;
  for( int i = 0; i < nrow; ++ i ) {
    energy += pow( yl[ i ], 2 );
  }
  energy /=  ( 2 * sigma2 );
  
  energy += log( sigma2 ) * nrow / (double)2;
  
  energy += ( a + 1 ) * log( mu - mk1 ) + b / ( mu - mk1 ) -
    log( pow( b, a ) / tgamma( a ) );
  
  energy += ( alphak + 1 ) * log( sigma2 ) + betak / sigma2 -
    log( pow( betak, alphak ) / tgamma( alphak ) );
  // Rprintf( "energy = %f\n", energy );
  double sum_theta = 0;
  if( nrow > 1) {
    for( int i = 0; i < 6; ++ i ) {
      sum_theta += pow( vtheta[ i ], 2 );
    }
  }
  
  energy += ( log( 2 * Pi ) * ncol / 2 + log( lambda2 ) * ncol / 2 +
    sum_theta / ( 2 * lambda2 ) ) / 2;
  // change ptr_seg[ 2, region ] back
  recoverLabel( region, ptr_seg );
  
  delete [] yln_;
  delete [] yln_i;
  delete [] yl_;
  delete [] yln;
  delete [] yl;
  delete [] vtheta;
  
  return energy;
}

// Single point energyY for healty voxel
double energyY( const int &curr_label,
                const int &curr_idx,
                const double &mu,
                const double &sigma2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta ) {
  double curr_intst = ptr_intst[ curr_idx - 1 ];
  curr_intst -= mu;
  double energy = curr_intst;
  double nbr_intst;
  
  for( int i = 0; i < 6; ++ i ) {
    int nidx =  ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nidx != NA_INTEGER ) {
      int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
      if( nlabel == curr_label ) {
        nbr_intst = ptr_nintst[ 6 * ( curr_idx - 1 ) + i ] - mu;
        energy -= theta[ i ] * nbr_intst;
      }
    }
  }
  energy = pow( energy, 2 ) / ( 2 * sigma2 ) + log( sigma2 ) / 2;
  return energy;
}

// Single point energyY for tumor voxel
// having tumor neighbor
double energyY( const int curr_label,
                const int curr_idx,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                const double &alphak,
                const double &betak,
                const double &a, const double &b ) {
  double curr_intst = ptr_intst[ curr_idx - 1 ];
  curr_intst -= mu;
  double energy = curr_intst;
  double nbr_intst;
  
  for( int i = 0; i < 6; ++ i ) {
    int nidx =  ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nidx != NA_INTEGER ) {
      int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
      if( nlabel == curr_label ) {
        nbr_intst = ptr_nintst[ 6 * ( curr_idx - 1 ) + i ] - mu;
        energy -= theta[ i ] * nbr_intst;
      }
    }
  }
  energy = pow( energy, 2 ) / ( 2 * sigma2 ) + log( sigma2 ) / 2;
  return energy;
}
// Single point energyY for tumor without tumor neighbor
// or outlier
double energyY( const int &curr_label,
                const int &curr_idx,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                const int *ptr_seg,
                const double *ptr_intst,
                const double &alphak,
                const double &betak,
                const double &a, const double &b ) {
  double curr_intst = ptr_intst[ curr_idx - 1 ];
  curr_intst -= mu;
  double energy = curr_intst;
  energy = pow( energy, 2 ) / ( 2 * sigma2 ) + log( sigma2 ) / 2;
  energy += ( a + 1 ) * log( mu - mk1 ) + b / ( mu - mk1 ) -
    log( pow( b, a ) / tgamma( a ) );
  energy += ( alphak + 1 ) * log( sigma2 ) + betak / sigma2 -
    log( pow( betak, alphak ) / tgamma( alphak ) );
  int ncol = 6;
  energy += ( log( 2 * Pi ) * ncol / 2 + log( lambda2 ) * ncol / 2 ) / 2;
  return energy;
}