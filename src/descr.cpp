#include <R.h>
#include <Rinternals.h>

#include <list>
#include <cmath>
#include <algorithm>

#include "descr.h"
#include "regions.h"
#include "calcDescr.h"

using std::list;

vector<double> descr2D( const int &len,
                        int *ptr_seg_copy, const int &label, 
                        vector<int> &region, 
                        const int *ptr_nidx, const int *ptr_aidx ) {
  vector<double> solidity( 3, 0 );
  vector<double> avg_spread( 3, 0 );
  vector<double> avg_round( 3, 0 );
  vector<double> res( 3, 0 ); // sodility, avg_spread, avg_round
  int size_region = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg_copy[ 2 * i ] == label ) {
      ++ size_region;
    }
  }
  // Find connected regions in ptr_seg
  vector<list<vector<int>>> slices = regions2D( ptr_seg_copy, len,
                                                region, label, 
                                                ptr_nidx,
                                                ptr_aidx );
  calcDescr( size_region,
             ptr_seg_copy, len, region, slices, 
             solidity, avg_spread, avg_round, 
             ptr_aidx, ptr_nidx );
  res[ 0 ] = *max_element( solidity.begin(), solidity.end() );
  res[ 1 ] = *max_element( avg_spread.begin(), avg_spread.end() );
  res[ 2 ] = *max_element( avg_round.begin(), avg_round.end() );
  Rprintf( "region size = %d; ", size_region );
  Rprintf( "sodility = %f, spread = %f, roundness = %f\n",
           res[ 0 ], res[ 1 ], res[ 2 ] );
  return res;
}
                           