#include "cmpE.h"
#include "label2col.h"
#include "getParm.h"

// compare energy for t1ce, flair or t2 images in estimation
void cmpE( const int &n_health,
            int idx, vector<double> &health_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma,
            vector<double> &theta ) {
  
  double energy;
  double min_energy;
  int min_label;
  double mu;
  double sigma2;
  int cidx;
  for( int i = - 1; i > - ( n_health + 1 ); -- i ) {
    cidx = label2col( i );
    getParm( mu, sigma2, theta, health_parm, cidx );
    energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                      ptr_intst, ptr_nintst, theta );
    energy += energyX( i, idx, false, ptr_seg, ptr_nidx,
                       ptr_delta, ptr_gamma[ 0 ] );
    if( i == - 1 ) {
      min_energy = energy;
      min_label = i;
    } else {
      if( energy < min_energy ) {
        min_energy = energy;
        min_label = i;
      }
    }
  }
  ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
}