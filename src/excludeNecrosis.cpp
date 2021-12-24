#include "excludeNecrosis.h"
#include "cnctRegion.h"
#include "zeroVector.h"
#include "pad2zero.h"

void excludeNecrosis( const int &len, const int *ptr_nidx,
                      int *ptr_new, 
                      const int *ptr_t1ce_csf,
                      const double *ptr_t1ce_intst,
                      int *ptr_add_necrosis,
                      vector<int> &add_necrosis_region,
                      const double &mean_csf ) {
  int idx = 0, n_necrosis = 0;
  double mean_necrosis = 0;
  zeroVector( ptr_add_necrosis, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_new[ 2 * i ] > 0 && 
        ptr_t1ce_csf[ 2 * i ] == 1 ) {
      ptr_add_necrosis[ 2 * i ] = 1;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_add_necrosis, 
                    ptr_add_necrosis, 1, add_necrosis_region ) ) {
      n_necrosis = add_necrosis_region.size();
      mean_necrosis = 0;
      for( vector<int>::const_iterator it = add_necrosis_region.begin();
           it != add_necrosis_region.end(); ++ it ) {
        idx = *it;
        mean_necrosis += ptr_t1ce_intst[ idx - 1 ]; 
      }
      mean_necrosis /= n_necrosis;
      if( mean_necrosis <= mean_csf ) {
        for( vector<int>::const_iterator it = add_necrosis_region.begin();
             it != add_necrosis_region.end(); ++ it ) {
          idx = *it;
          ptr_new[ 2 * ( idx - 1 ) ] = 0;
        }
      }
    } 
  }
  pad2zero( ptr_add_necrosis, len );
}