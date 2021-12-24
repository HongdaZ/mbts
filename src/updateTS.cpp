#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <list>
#include <vector>

#include "updateTS.h"
#include "initMV.h"

using std::list;
using std::vector;

// region starts from 1
// updateTheta and sigma2 for health and tumorous regions
void updateTS( const int &nrow,
               const double &mu,
               double &sigma2,
               const double &lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               vector<double> &theta,
               const double &alphal,
               const double &betal,
               double *yl, double *yln, double *ylna,
               const double *yln_,
               const double *yln_i,
               const double *yl_ ){
  int ncol = 6 / 2; // only 6 / 2 different theta's
  // Initialize yln and yl;
  initMV( yl, yl_, yln, ylna, yln_, yln_i, nrow, mu );
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "%f\n", yl[ i ] );
  // }
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "\n" );
  //   for( int j = 0; j < 6; ++ j ) {
  //     Rprintf( "%f,\t", yln[ nrow * j + i ] );
  //   }
  // }
  
  // Identity matrix
  double *I = new double[ ncol * ncol ];
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = 0; j < ncol; ++ j ) {
      I[ i * ncol + j ] = ( i == j ) ? 1 : 0;
    }
  }
  double *tmp = new double[ ncol ];
  double *vthetah = new double[ ncol ]; 
  // t( Ylna ) Ylna / sigma2 + I / lambda2 
  double alpha = 1 / sigma2;
  double beta = 1 / lambda2;
  // Rprintf( "nrow = %d, ncol = %d, alpha = %f, beta = %f\n", 
  //          nrow, ncol, alpha, beta );
  F77_CALL( dsyrk )( "u", "t", &ncol, &nrow, &alpha, ylna, &nrow, &beta, 
            I, &ncol );
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = i + 1; j < ncol; ++ j ) {
      I[ i * ncol + j ] = I[ j * ncol + i ];
    }
  }
  
  // inverse
  int *ipiv = new int[ ncol ];
  int lwork = ncol * ncol;
  double *work = new double[ lwork ];
  int info;
  F77_CALL( dgetrf )( &ncol, &ncol, I, &ncol, ipiv, &info );
  F77_CALL( dgetri )( &ncol, I, &ncol, ipiv, work, &lwork, &info );
  // // Debug for inverse 
  // for( int i = 0; i < ncol; ++ i ) {
  //   for( int j = 0; j < ncol; ++ j ) {
  //     Rprintf( "%f\t", I[ i * ncol + j ] );
  //   }
  //   Rprintf( "\n" );
  // }
  // t( ylna ) yl \ sigma2
  int incx = 1;
  double *x = new double[ nrow ];
  beta = 0;
  F77_CALL( dgemv )( "t", &nrow, &ncol, &alpha, ylna, &nrow, yl, &incx, &beta, 
            tmp, &incx );
  // // debug matrix vector multiply
  // for( int i = 0; i < ncol; ++ i ) {
  //   Rprintf( "%f\t", tmp[ i ] );
  // }
  // ////////////////////////////////////
  alpha = 1;
  F77_CALL( dgemv )( "t", &ncol, &ncol, &alpha, I, &ncol, tmp, &incx, &beta, 
            vthetah, &incx );
  for( int i = 0; i < ncol; ++ i ) {
    theta[ i ] = vthetah[ i ];
    theta[ 2 * ncol - 1 - i ] = vthetah[ i ];
  }
  // // debug matrix vector multiply
  // for( int i = 0; i < ncol; ++ i ) {
  //   Rprintf( "%f\t", theta[ i ] );
  // }
  // ////////////////////////////////////
  
  // update sigma2
  // yl - yln * theta; = yl - ylna * thetah
  alpha = - 1;
  beta = 1;
  F77_CALL( dgemv )( "n", &nrow, &ncol, &alpha, ylna, &nrow, vthetah, &incx, 
            &beta, yl, &incx );
  // // debug matrix vector multiply
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "%f\t", yl[ i ] );
  // }
  // Rprintf( "\n");
  // ////////////////////////////////////
  double sum = betal;
  for( int i = 0; i < nrow; ++ i ) {
    sum += pow( yl[ i ], 2. ) / 2.;
  }
  sigma2 = sum / ( nrow / 2. + alphal + 1 );

  delete [] I;
  delete [] tmp;
  delete [] vthetah;
  delete [] ipiv;
  delete [] work;
  delete [] x;
  
  
  return;
  
}
