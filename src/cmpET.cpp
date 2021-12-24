#include <R.h>
#include <Rinternals.h>

#include "cmpET.h"
#include "findOutLabel.h"
#include "newTumorLabel.h"
#include "eraseRegion.h"
#include "eraseOutl.h"
#include "addRegion.h"
#include "addOutl.h"
#include "tumorNbr.h"
#include "getRegion.h"
#include "assignParm.h"
#include "label2col.h"
#include "getParm.h"
#include "zeroVector.h"
#include "addVoxel.h"
#include "eraseVoxel.h"
#include "clearVector.h"

// compare energy for training
void cmpET( const int &n_health,
            const list<int> &update_parm,
            vector<int> &region, const int &idx, const int &sc,
            const vector<int> &regions_whole,
            const vector<int> &regions_sub,
            vector<int> &tumor_labels, vector<int> &outl_labels,
            const vector<double> &health_parm,
            vector<double> &tumor_parm, 
            vector<double> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2,
            // outlier_parm( 3, 0 )
            vector<double> &outlier_parm,
            // theta( 6, 0 )
            vector<double> &theta,
            // tmp_parm( 9, 0 )
            vector<double> &tmp_parm,
            // out_theta( 6, 0 )
            vector<double> &out_theta,
            // new_out_parm( 2, 0 )
            vector<double> &new_out_parm,
            // whole_parm( 8, 0 )
            vector<double> &whole_parm,
            // label_whole_parm( 9, 0 )
            vector<double> &label_whole_parm,
            // used to store parameters of whole region and sub-regions
            vector<double> &region_parm,
            int &n_tumor, int &n_outl,
            const int &n_region,
            list<list<int>> &tumor_regions,
            const int &n_row,
            const list<int> &tumor_label ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  int start = -1;
  double mu, sigma2, energy;
  int cidx = 0;
  if( sc != 0 ) {  // split or combine
    // energy of whole region, subregions and outlier
    vector<double> nrg( n_region + 1, 0 );
    // parameters for regions
    clearVector( region_parm );
    region_parm.resize( n_row * n_region );
    
    int region_label;
    int region_idx;
    sigma2 = ptr_beta[ n_health ] / ( ptr_alpha[ n_health ] + 1 );
    list<int>::const_iterator it = update_parm.begin();
    for( int i = 0; i < n_region; ++ i ) {
      // get region and region label
      getRegion( region_label, region, regions_whole, regions_sub,
                 start );
      if( region.size() == 1 ) {
        region_idx = region[ 0 ];
        if( *it == 0 ) {
          updateParm( mu, sigma2, region_idx, ptr_m[ n_health ], 
                      ptr_m[ n_health - 1 ],
                      ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, ptr_seg,
                      ptr_alpha[ n_health ], ptr_beta[ n_health ], 20 );
        } else if( *it < - 3 ){
          cidx = label2col( *it );
          getParm( mu, sigma2, theta, tumor_parm, cidx );
        }
        // region_label, mu, sigma2, theta
        region_parm[ n_row * i ] = region_label;
        region_parm[ n_row * i + 1 ] = mu;
        region_parm[ n_row * i + 2 ] = sigma2;
        for( int j = 0; j < 6; ++ j ) {
          region_parm[ n_row * i + j + 3 ] = 0;
        }
        // calculate energy
        energy = energyY( region_label, region_idx, mu,
                          ptr_m[ n_health - 1 ], sigma2,
                          ptr_lambda2[ n_health ], ptr_seg,
                          ptr_intst, ptr_alpha[ n_health ], 
                          ptr_beta[ n_health ], ptr_a[ 0 ], ptr_b[ 0 ] );
      } else {
        if( *it == 0 ) {
          updateParm( mu, theta, sigma2, region, ptr_m[ n_health ], 
                      ptr_m[ n_health - 1 ], ptr_a[ 0 ], ptr_b[ 0 ],
                      ptr_intst, region_label, ptr_lambda2[ n_health ],
                      ptr_seg, ptr_nidx, ptr_nintst, ptr_alpha[ n_health ],
                      ptr_beta[ n_health ], 1 );
        } else if( *it < -3 ){
          cidx = label2col( *it );
          getParm( mu, sigma2, theta, tumor_parm, cidx );
        }
        // region_label, mu, sigma2, theta
        region_parm[ n_row * i ] = region_label;
        region_parm[ n_row * i + 1 ] = mu;
        region_parm[ n_row * i + 2 ] = sigma2;
        for( int j = 0; j < 6; ++ j ) {
          region_parm[ n_row * i + j + 3 ] = theta[ j ];
        }
        // calculate energy
        energy = energyY( region, mu, ptr_m[ n_health - 1 ], sigma2,
                          ptr_lambda2[ n_health ], ptr_seg, ptr_nidx,
                          ptr_intst, ptr_nintst,
                          theta, ptr_alpha[ n_health ], 
                          ptr_beta[ n_health ],
                          ptr_a[ 0 ], ptr_b[ 0 ] );
      }
      nrg[ i ] = energy;
    }
    
    // New outlier label
    int out_label;
    findOutLabel( out_label, idx, ptr_seg, outl_labels );
    
    // outlier parameters
    double &out_mu = mu;
    double &out_sigma2 = sigma2;
    double &out_energy = energy;
    out_sigma2 = ptr_beta[ n_health ] / ( ptr_alpha[ n_health ] + 1 );
    updateParm( out_mu, out_sigma2, idx, ptr_m[ n_health ], 
                ptr_m[ n_health - 1 ], ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst,
                ptr_seg, ptr_alpha[ n_health ], ptr_beta[ n_health ], 20 );
    outlier_parm[ 0 ] = out_label;
    outlier_parm[ 1 ] = out_mu;
    outlier_parm[ 2 ] = out_sigma2;
    
    out_energy = energyX( out_label, idx, true, ptr_seg, ptr_nidx,
                          ptr_delta, ptr_gamma[ 0 ] );
    out_energy += energyY( out_label, idx, out_mu, ptr_m[ n_health - 1 ],
                           out_sigma2,
                           ptr_lambda2[ n_health ], ptr_seg, ptr_intst,
                           ptr_alpha[ n_health ],
                           ptr_beta[ n_health ], ptr_a[ 0 ], ptr_b[ 0 ] );
    nrg.back() = out_energy;
    double combine_nrg = nrg[ 0 ];
    double split_nrg = 0;
    for( int i = 1; i < ( n_region + 1 ); ++ i ) {
      split_nrg += nrg[ i ];
    }
    for( int i = 0; i < n_row; ++ i ) {
      label_whole_parm[ i ] = region_parm[ i ];
    }
    int whole_label = label_whole_parm[ 0 ];
    // update parameters for whole region
    if( combine_nrg <= split_nrg && sc == 1 ) {
      for( int i = 0; i < 8; ++ i ) {
        whole_parm[ i ] = label_whole_parm[ 1 + i ];
      }
      assignParm( tumor_parm, whole_label, whole_parm );
    // remove old whole region and add new splitted regions and new outlier
    } else if ( combine_nrg > split_nrg && sc == 1 ) {
      eraseRegion( whole_label, tumor_labels, tumor_parm,
                   tumor_regions, n_tumor );
      vector<double> new_region_parm( region_parm.begin() + n_row,
                                      region_parm.end() );
      // add sub-regions
      addRegion( ptr_seg, new_region_parm, n_row, tumor_labels,
                 tumor_parm, tumor_regions, regions_sub, n_tumor );
      for( int i = 0; i < 2; ++ i ) {
        new_out_parm[ i ] = outlier_parm[ 1 + i ];
      }
      addOutl( ptr_seg, idx, out_label, new_out_parm,
               outl_labels, outl_parm,
               n_outl );
      // if( out_label == 32 ) {
      //   Rprintf( "place 1: idx = %d; mu = %f, sigma2 = %f\n", idx, new_out_parm[ 0 ], 
      //            new_out_parm[ 1 ] );
      // }
      // remove old subregions, parameters and labels,
      // add new whole regions label,region and parameters,
      // set seg to tumor,
      // remove outlier from outl_labels
      // remove outl_parm
    } else if( combine_nrg < split_nrg && sc == 2 ) {
      for( list<int>::const_iterator it = tumor_label.begin();
           it != tumor_label.end(); ++ it ) {
        eraseRegion( *it, tumor_labels, tumor_parm, tumor_regions,
                     n_tumor );
      }
      // add whole region
      addRegion( ptr_seg, label_whole_parm, n_row, tumor_labels,
                 tumor_parm, tumor_regions, regions_whole, n_tumor );
      eraseOutl( out_label, outl_labels, outl_parm, n_outl );
      // if( out_label == 32 ) {
      //   bool in = outl_parm.find( 32 ) != outl_parm.end();
      //   Rprintf( "place 5: idx = %d, in list = %d\n", 
      //            idx, in );
      // }
      
    // update parameters for subregions
    } else if( combine_nrg >= split_nrg && sc == 2 ) {
      int sub_label;
      for( int i = 1; i < n_region; ++ i ) {
        sub_label = region_parm[ n_row * i ];
        vector<double> &new_parm = whole_parm;
        for( int j = 0; j < 8; ++ j ) {
          new_parm[ j ] = region_parm[ n_row * i + j + 1 ];
        }
        assignParm( tumor_parm, sub_label, new_parm );
      }
    }
    // no split or combine
  } else {
    double min_energy;
    int min_label;
    if( curr_label > - 4 && curr_label < 1 ) {
      double &h_energy = energy;
      double &h_mu = mu, &h_sigma2 = sigma2;
      vector<double> &h_theta = theta;
      for( int i = - 1; i > - ( n_health + 1 ); -- i ) {
        cidx = label2col( i );
        getParm( h_mu, h_sigma2, h_theta, health_parm, cidx );
        h_energy = energyY( i, idx, h_mu, h_sigma2, ptr_seg, ptr_nidx,
                            ptr_intst, ptr_nintst, h_theta );
        h_energy += energyX( i, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta, ptr_gamma[ 0 ] );
        // Rprintf( "h_label = %d; h_energy = %f; h_mu = %f; h_sigma2 = %f;\n",
        //          i, energy, mu, sigma2 );
        if( i == - 1 ) {
          min_energy = h_energy;
          min_label = i;
        } else {
          if( h_energy < min_energy ) {
            min_energy = h_energy;
            min_label = i;
          }
        } 
      }
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
      // tumor or outlier
    } else {
      vector<double> &t_theta = theta;
      double t_energy;
      // check if having tumor neighbor
      int tumor_idx = firstTumorNbr( idx, ptr_seg, ptr_nidx ) ;
      int t_label;
      // have tumor neighbor
      if( tumor_idx != 0 ) {
        t_label = ptr_seg[ 2 * ( tumor_idx - 1 ) ];
        // tumor energy
        cidx = label2col( t_label );
        double &t_mu = mu, &t_sigma2 = sigma2;
        getParm( t_mu, t_sigma2, t_theta, tumor_parm, cidx );
        t_energy = energyY( t_label, idx, t_mu, ptr_m[ n_health - 1 ],
                            t_sigma2, ptr_lambda2[ n_health ], ptr_seg,
                            ptr_nidx, ptr_intst, ptr_nintst,
                            t_theta, ptr_alpha[ n_health ], 
                            ptr_beta[ n_health ], ptr_a[ 0 ], ptr_b[ 0 ] );
        t_energy += energyX( t_label, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta, ptr_gamma[ 0 ] );
        // Rprintf( "t_label = %d; t_energy = %f; t_mu = %f; t_sigma2 = %f;\n",
        //          t_label, t_energy, t_mu, t_sigma2 );
        // Outlier energy
        // New outlier label
        int out_label;
        findOutLabel( out_label, idx, ptr_seg, outl_labels );
        double &out_energy = energy;
        out_energy = energyX( out_label, idx, false, ptr_seg,
                              ptr_nidx, ptr_delta,
                              ptr_gamma[ 0 ] );
        // outlier parameters
        double &out_mu = mu, &out_sigma2 = sigma2;
        if( curr_label > 0 ) {
          cidx = label2col( curr_label );
          getParm( out_mu, out_sigma2, outl_parm, cidx );
        } else {
          out_sigma2 = ptr_beta[ n_health ] / ( ptr_alpha[ n_health ] + 1 );
          updateParm( out_mu, out_sigma2, idx, ptr_m[ n_health ],
                      ptr_m[ n_health - 1 ], ptr_a[ 0 ], ptr_b[ 0 ],
                      ptr_intst, ptr_seg, ptr_alpha[ n_health ], 
                      ptr_beta[ n_health ], 20 );
          new_out_parm[ 0 ] = out_mu;
          new_out_parm[ 1 ] = out_sigma2;
        }
        out_energy += energyY( out_label, idx, out_mu, 
                               ptr_m[ n_health - 1 ],
                               out_sigma2, ptr_lambda2[ n_health ], 
                               ptr_seg, ptr_intst, ptr_alpha[ n_health ],
                               ptr_beta[ n_health ],
                               ptr_a[ 0 ], ptr_b[ 0 ] );
        // Rprintf( "out_label = %d; out_energy = %f; out_mu = %f; out_sigma2 = %f\n",
        //          out_label, out_energy, out_mu, out_sigma2 );
        if( curr_label > 0 && t_energy <= out_energy ) {
          addVoxel( idx, min_label, tumor_regions, ptr_seg );
          eraseOutl( out_label, outl_labels, outl_parm, n_outl );
          // if( curr_label == 32 ) {
          //   bool in = outl_parm.find( 32 ) != outl_parm.end();
          //   Rprintf( "place 3: idx = %d, in list = %d\n", idx, in );
          // }
        } else if( curr_label < - 3 && t_energy >= out_energy ) {
          eraseVoxel( idx, curr_label, tumor_regions );
          addOutl( ptr_seg, idx, min_label, new_out_parm,
                   outl_labels, outl_parm, n_outl );
          // if( out_label == 32 ) {
          //   Rprintf( "place 2: idx = %d; mu = %f, sigma2 = %f\n", idx, new_out_parm[ 0 ], 
          //            new_out_parm[ 1 ] );
          // }
        }
        // not have tumor  
      } else {
        if( curr_label > 0 ) {
          // Find a new tumor region label;
          int new_label;
          int start = 0;
          newTumorLabel( new_label, start, tumor_labels );
          // new tumor region parameters
          double &t_mu = mu, &t_sigma2 = sigma2;
          t_sigma2 = ptr_beta[ n_health ] / ( ptr_alpha[ n_health ] + 1 );
          // use the function for outliers and single voxel tumor regions
          updateParm( t_mu, t_sigma2, idx, ptr_m[ n_health ],
                      ptr_m[ n_health - 1 ], ptr_a[ 0 ], ptr_b[ 0 ],
                      ptr_intst, ptr_seg,
                      ptr_alpha[ n_health ],
                      ptr_beta[ n_health ], 20 );
          vector<double> &new_parm = tmp_parm;
          new_parm[ 0 ] = t_label;
          new_parm[ 1 ] = t_mu;
          new_parm[ 2 ] = t_sigma2;
          for( int i = 0; i < 6; ++ i ) {
            new_parm[ i + 3 ] = 0;
          }
          eraseOutl( curr_label, outl_labels, outl_parm, n_outl );
          addRegion( ptr_seg, new_parm, idx, tumor_labels,
                     tumor_parm, tumor_regions, n_tumor );
          // if( curr_label == 32 ) {
          //   bool in = outl_parm.find( 32 ) != outl_parm.end();
          //   Rprintf( "place 4: idx = %d, in list = %d\n", idx, in );
          // }
        }
      }
    }
  }
}