#include <R.h>
#include <Rinternals.h>

#include "initParm.h"
#include "updateParm.h"
#include "getRegion.h"
#include "clearVector.h"
#include "assignParm.h"
#include "label2col.h"
#include "zeroVector.h"

// Initialize parameters
void initParm( const int &n_health,
               vector<int> &region, vector<double> &theta, 
               const bool &first_run, vector<double> &health_parm,
               vector<double> &tumor_parm,
               int *ptr_seg, const double *ptr_m,
               const double *ptr_nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               const list<list<int>> &tumor_regions,
               const double *ptr_a, const double *ptr_b, const int &len, 
               const int &maxit ) {
  // update parameters for tumor regions
  double mu, sigma2;
  int curr_label;
  list<list<int>>::const_iterator it = tumor_regions.begin();
  for( ;it != tumor_regions.end(); ++ it ) {
    curr_label = it->front();
    sigma2 = ptr_beta[ n_health ] / ( ptr_alpha[ n_health ] + 1 );
    // region starts from 1
    const list<int> &t_region = *it;
    getRegion( region, t_region );
    if( first_run ) {
      if( region.size() == 1 ) {
        zeroVector( theta );
        updateParm( mu, sigma2, region.front(), ptr_m[ n_health ], 
                    ptr_m[ n_health - 1 ], 
                    ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst,
                    ptr_seg, ptr_alpha[ n_health ],
                    ptr_beta[ n_health ], maxit );
      } else {
        updateParm( mu, theta, sigma2, region, ptr_m[ n_health ], 
                    ptr_m[ n_health - 1 ],
                    ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst,
                    curr_label, ptr_lambda2[ n_health ],
                    ptr_seg, ptr_nidx, 
                    ptr_nintst, ptr_alpha[ n_health ],
                    ptr_beta[ n_health ], maxit );
      }
      assignParm( tumor_parm, curr_label, mu, sigma2, theta );
    } else if( it->size() > 27 ) {
      updateParm( mu, theta, sigma2, region, ptr_m[ n_health ], 
                  ptr_m[ n_health - 1 ], ptr_a[ 0 ], ptr_b[ 0 ], 
                  ptr_intst, curr_label, ptr_lambda2[ n_health ], 
                  ptr_seg, ptr_nidx, ptr_nintst, ptr_alpha[ n_health ], 
                  ptr_beta[ n_health ], maxit );
      assignParm( tumor_parm, curr_label, mu, sigma2, theta );
    }
  }
  // Rprintf( "update parm for healthy regions\n" );
  // Initialize parameters for healthy regions
  int h_idx;
  for( int i = - 1; i > - ( n_health + 1 );  -- i ) {
    curr_label = i;
    h_idx = label2col( curr_label );// == 0, 1, 2
    sigma2 = ptr_beta[ h_idx ] / ( ptr_alpha[ h_idx ] + 1 );
    // Rprintf( "curr_label = %d \n", curr_label );
    getRegion( region, curr_label, ptr_seg, len );
    // Rprintf( "curr_label = %d, region size = %d\n", 
    //          curr_label, region.size() );
    if( region.size() == 0 ) {
      continue;
    }
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], 
                ptr_nu2[ h_idx ], ptr_intst, curr_label, 
                ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx, ptr_nintst,
                ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );
    assignParm( health_parm, curr_label, mu, sigma2, theta );
  }
  return;
}
// Initialize parameters for t1ce, flair and t2 images
void initParmHealth( const int &n_health, 
                     vector<int> &region, vector<double> &theta,
                     vector<double> &health_parm,
                     int *ptr_seg, const double *ptr_m,
                     const double *ptr_nu2, const double *ptr_intst, 
                     const double *ptr_lambda2, const int *ptr_nidx,
                     const double *ptr_nintst, const double *ptr_alpha,
                     const double *ptr_beta, const int &len, 
                     const int &maxit ) {
  // Initialize parameters for healthy regions
  int curr_label;
  int h_idx;
  double mu, sigma2;
  for( int i = - 1; i > - ( n_health + 1 );  -- i ) {
    curr_label = i;
    h_idx = label2col( curr_label );// == 0, 1, 2 or 0, 1
    sigma2 = ptr_beta[ h_idx ] / ( ptr_alpha[ h_idx ] + 1 );
    getRegion( region, curr_label, ptr_seg, len );
    if( region.size() == 0 ) {
      continue;
    }
    h_idx = label2col( curr_label ); // == 0, 1, 2
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], 
                ptr_nu2[ h_idx ], ptr_intst, curr_label,
                ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx, ptr_nintst, 
                ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );
    
    assignParm( health_parm, curr_label, mu, sigma2, theta );
  }
  return;
}