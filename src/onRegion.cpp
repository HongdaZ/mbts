#include <R.h>
#include <Rinternals.h>

#include <cmath>
#include <algorithm>

#include "onRegion.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "zeroVector.h"
#include "pad2zero.h"
#include "spread.h"
#include "pTNbr.h"
#include "inRegion.h"
#include "clearVector.h"
#include "radius.h"
#include "expRatio.h"
#include "areaRatio.h"
#include "ballCrownVol.h"
#include "thick.h"

using std::cbrt;
using std::max_element;

// Extend ptr_seg1 to adjacent regions of ptr_seg2
void onRegion( int *ptr_on, const int &len, const double &prop,
               int *ptr_seg1, const int &label1, 
               int *ptr_seg2, const int &label2,
               int *ptr_seg2_copy,
               vector<int> &region, const double &spread_factor,
               const int *ptr_nidx, const int *ptr_aidx,
               const int &nr, const int &nc, const int &ns,
               int *ptr_seg_copy, const double &prop_hull,
               const double &nt_prop ) {
  vector<double> volume;
  double max_volume;
  double prop_tumor;
  int idx = 0, n_idx = 0;
  double p_t_nbr = 0;
  int n_enclose_t;
  double spread_idx = 0;
  bool add = false;
  double r = 0;
  int n_cnct_old;
  int n_new = 0;
  double h, v, area_ratio, exp_ratio, r0, roundness;
  
  int *ptr_enclose_tumor = new int[ 2 * len ]();
  int *ptr_new_region = new int[ 2 * len ]();
  int *ptr_tmp_seg1 = new int[ 2 * len ]();
  int *ptr_old_new = new int[ 2 * len ]();
  int *ptr_cnct_old = new int[ 2 * len ]();
  
  for( int i = 0; i < len; ++ i ) {
    ptr_tmp_seg1[ 2 * i ] = ptr_seg1[ 2 * i ]; 
  }
  
  vector<int> tmp_region;
  tmp_region.reserve( len );
  
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_seg2, ptr_seg2,
                    label2, region ) ) {
      int n_idx = excldRegion( region, ptr_nidx, ptr_seg2,
                               ptr_seg1, label1 );
      if(  n_idx != 0  ) {
        add = false;
        zeroVector( ptr_new_region, len );
        zeroVector( ptr_enclose_tumor, len );
        zeroVector( ptr_old_new, len );
        zeroVector( ptr_cnct_old, len );
        n_new = region.size();
        
        for( vector<int>::const_iterator it = region.begin();
             it != region.end(); ++ it ) {
          idx = *it;
          if( idx != 0 ) {
            ptr_new_region[ 2 * ( idx - 1 ) ] = 1;
          }
        }
        for( int j = 0; j < len; ++ j ) {
          if( ptr_new_region[ 2 * j ] == 1 ||
              ptr_tmp_seg1[ 2 * j ] == label1 ) {
            ptr_old_new[ 2 * j ] = 1;
          }
        }
        cnctRegion( i + 1, ptr_nidx, ptr_old_new, ptr_old_new,
                    1, tmp_region );
        zeroVector( ptr_old_new, len );
        n_cnct_old = 0;
        for( vector<int>::const_iterator it_tmp = tmp_region.begin();
             it_tmp != tmp_region.end(); ++ it_tmp ) {
          idx = *it_tmp;
          ptr_old_new[ 2 * ( idx - 1 ) ] = 1;
          if( ptr_new_region[ 2 * ( idx - 1 ) ] != 1 ) {
            ptr_cnct_old[ 2 * ( idx - 1 ) ] = 1;
            ++ n_cnct_old;
          } else {
            ptr_cnct_old[ 2 * ( idx - 1 ) ] = 0;
          }
        }
        volume = inRegion( ptr_enclose_tumor, len, ptr_new_region, 1,
                           ptr_cnct_old, 1, 
                           ptr_seg2_copy,
                           tmp_region, ptr_nidx,
                           ptr_aidx, nr, nc, ns );
        max_volume = *max_element( volume.begin(), volume.end() );
        n_enclose_t = 0;
        for( int j = 0; j < len; ++ j ) {
          if( ptr_enclose_tumor[ 2 * j ] == 1 ) {
            ++ n_enclose_t;
          }
        }
        prop_tumor = ( double ) ( n_new + n_enclose_t ) / max_volume;
        // Rprintf( "max_volume = %f, ", max_volume );
        // Rprintf( "prop_tumor = %f, prop_hull = %f\n",
        //          prop_tumor, prop_hull );
        if( prop_tumor > prop_hull ) {
          // Rprintf( "n_enclose_t = %d, nt_prop * n_new = %f\n",
          //          n_enclose_t, nt_prop * n_new );
          if( n_enclose_t > nt_prop * n_new ) {
            clearVector( tmp_region );
            for( int j = 0; j < len; ++ j ) {
              if( ptr_enclose_tumor[ 2 * j ] == 1 ||
                  ptr_new_region[ 2 * j ] == 1 ) {
                tmp_region.push_back( j + 1 );
              }
            }
            spread_idx = spread( tmp_region, ptr_seg_copy, 
                                 len, ptr_nidx );
            // Rprintf( "spread_idx = %f, spread_add = %f\n",
            //          spread_idx, spread_factor );
            if( spread_idx < spread_factor ) {
              add = true;
            } else {
              add = false;
            }
          } else {
            p_t_nbr = pTNbr( region, ptr_seg1, label1, ptr_nidx );
            r0 = radius2D( p_t_nbr * region.size() ); 
            h = thick( ptr_tmp_seg1, label1, ptr_new_region, 1,
                       len, ptr_nidx );
            area_ratio = areaRatio( ptr_tmp_seg1, label1, 
                                    ptr_new_region, 1, 
                                    len, ptr_nidx );
            
            exp_ratio = expRatio( r0, h );
            v = ballCrownVol( r0, h );
            if( region.size() > v * 0.75 &&
                area_ratio > exp_ratio * 0.75 ) {
              spread_idx = spread( region, ptr_seg_copy, 
                                   len, ptr_nidx );
              
              
              if( spread_idx < spread_factor ) {
                add = true;
              } else {
                add = false;
              }
            } else {
              add = false;
            }
          }
        } else {
          add = false;
        }
        
        if( add ) {
          if( n_new < prop * n_cnct_old ) {
            for( vector<int>::const_iterator it = region.begin();
                 it != region.end(); ++ it ) {
              idx = *it;
              if( idx != 0 ) {
                ptr_on[ 2 * ( idx - 1 ) ] = 1;
              }
            }
          }
        }
      }
    }
  }
  
  delete [] ptr_enclose_tumor;
  delete [] ptr_new_region;
  delete [] ptr_tmp_seg1;
  delete [] ptr_old_new;
  delete [] ptr_cnct_old;
  
  pad2zero( ptr_seg2, len );
}