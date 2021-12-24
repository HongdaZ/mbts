#include <R.h>
#include <Rinternals.h>



#include "scTrn.h"
#include "findRegion.h"
#include "newTumorLabel.h"
#include "tumorNbr.h"
#include "nonTumor.h"
#include "clearVector.h"
#include "zeroVector.h"
#include "assignRegion.h"
#include "getRegion.h"

using std::list;
using std::vector;

// start = 1 to length array
// nidx also starts from 1
int scTrn( int &n_region, list<int> &update_parm,
           list<list<int>> &tumor_regions,
           vector<int> &region,
           const vector<int> &tumor_labels,
           int *ptr_label, const int *ptr_nidx,
           const int &len, const int &start,
           vector<int> &regions_whole, vector<int> &regions_sub,
           list<int> &tumor_nbr, list<int> &tumor_label ) {
  list<double> region_size;
  // 0: updateParm; < 0; getParm
  clearVector( update_parm );
  
  int current = ptr_label[ 2 * ( start - 1 ) ];
  if( current <= 0 && current >= -3 ) {
    return 0;
  } else {
    n_region = 0;
    const double proportion = .95;
    clearVector( regions_whole );
    clearVector( regions_sub );
    clearVector( tumor_nbr );
    clearVector( tumor_label );
    int current = ptr_label[ 2 * ( start - 1 ) ];
    
    // Add tumor neighbors to list
    tumorNbrLabel( tumor_label, tumor_nbr, start, ptr_label, ptr_nidx );
    
    if( tumor_nbr.size() < 2 ) {
      return 0;
      // possibly split
    } else if( current <= - 4 ) {
      regions_whole.push_back( current );
      regions_whole.push_back( start );
      // remember to reset seg[ 2, ] == 0
      ptr_label[ 2 * start - 1  ] = 1;
      // remember to reset seg[ 2, ] == 0
      for( list<int>::iterator it = tumor_nbr.begin();
           it != tumor_nbr.end(); ++ it ) {
        ptr_label[ 2 * *it - 1 ] = 2;
      }
      bool early_return = 0;
      int new_label_start = 0;
      int new_start;
      // find sub-regions
      while( !tumor_nbr.empty() ) {
        new_start = tumor_nbr.front();
        tumor_nbr.pop_front();
        ptr_label[ 2 * new_start - 1 ] = 0;
        findRegion( region, n_region, ptr_label, ptr_nidx,
                    false, tumor_nbr, early_return, new_start );
        if( early_return ) {
          ptr_label[ 2 * start - 1  ] = 0;
          return 0;
        }
        int region_label;
        if( n_region == 0 ) {
          region_label = current;
        } else {
          newTumorLabel( region_label, new_label_start, tumor_labels );
          
        }
        assignRegion( regions_whole, regions_sub, region, region_label );
        region_size.push_back( region.size() );
        ++ n_region;
      }
      region_size.push_front( regions_whole.size() - 1 );
      regions_whole.push_back( NA_INTEGER );
      // restore the value of seg[ 2, ] to 0;
      ptr_label[ 2 * start - 1  ] = 0;
      if( n_region < 2 ) {
        return 0;
      } else {
        // count whole region
        ++ n_region;
        list<double>::iterator it = region_size.begin();
        int whole_size = *it;
        // getParm or updateParm
        if( whole_size > 100 ) {
          update_parm.push_back( current );
          for( ++ it; it != region_size.end(); ++ it ) {
            if( ( *it / whole_size ) > proportion ) {
              update_parm.push_back( current );
            } else {
              update_parm.push_back( 0 );
            }
          }
        } else {
          for( int i = 0; i < n_region; ++ i ) {
            update_parm.push_back( 0 );
          }
        }
        return 1;
      }
    } else {
      if( tumor_label.size() == 1 ) {
        return 0;
      }
      // possibly combine
      int combine_label = tumor_label.back();
      regions_whole.push_back( combine_label );
      regions_whole.push_back( start );
      // add whole region
      ++ n_region;
      // add sub-regions
      n_region += tumor_label.size();
      int r_size;
      for( list<int>::iterator it = tumor_label.begin();
           it != tumor_label.end(); ++ it ) {
        assignRegion( regions_whole, regions_sub, tumor_regions,
                      *it, r_size );
        region_size.push_back( r_size );
      }
      region_size.push_front( regions_whole.size() - 1 );
      regions_whole.push_back( NA_INTEGER );
      
      list<double>::iterator it = region_size.begin();
      int whole_size = *it;
      if( whole_size > 100 ) {
        int combine_parm_label = 0;
        update_parm = tumor_label;
        list<int>::const_iterator it_tlabel = tumor_label.begin();
        for( ++ it; it != region_size.end(); ++ it, ++ it_tlabel ) {
          if( ( *it / whole_size ) > proportion ) {
            combine_parm_label = *it_tlabel;
            break;
          }
        }
        update_parm.push_front( combine_parm_label );
      } else {
        for( int i = 0; i < n_region; ++ i ) {
          update_parm.push_back( 0 );
        }
      }
      return 2;
    }
  }
}