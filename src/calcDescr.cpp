#include <R.h>
#include <Rinternals.h>

#include "calcDescr.h"
#include "cHull.h"
#include "polygonArea.h"
#include "spread.h"
#include "roundness.h"
#include "clearVector.h"

void calcDescr( const int &size_region,
                int *ptr_seg_copy, const int &len, 
                vector<int> &region,
                const vector<list<vector<int>>> &slices, 
                vector<double> &solidity, 
                vector<double> &avg_spread, 
                vector<double> &avg_round,
                const int *ptr_aidx, const int *ptr_nidx ) {
  // Sagittal, coronal and axial
  vector<int> row{ 2, 3, 1, 3, 1, 2 };
  vector<double> shift{ -0.5, 0.5, 
                       -0.5, -0.5, 
                       0.5, -0.5, 
                       0.5, 0.5 };
  double hull_area, c_spread, c_roundness;
  int len_slice;
  double c_1, c_2;
  for( int i = 0; i < 3; ++ i ) {
    solidity[ i ] = 0;
    avg_spread[ i ] = 0;
    avg_round[ i ] = 0;
    const list<vector<int>> &p_slice = slices[ i ];
    for( list<vector<int>>::const_iterator it_s = p_slice.begin();
         it_s != p_slice.end(); ++ it_s ) {
      
      len_slice = it_s->size() / 4;
      vector<double> chull_input( 2 * len_slice * 4 );
      clearVector( region );
      for( int j = 0; j < len_slice; ++ j ) {
        region.push_back( ( *it_s )[ 4 * j ] );
        c_1 = ( *it_s )[ 4 * j + row[ 2 * i ] ];
        c_2 = ( *it_s )[ 4 * j + row[ 2 * i + 1 ] ];
        for( int k = 0; k < 4; ++ k ) {
          chull_input[ 8 * j + 2 * k ] = c_1 + shift[ 2 * k ];
          chull_input[ 8 * j + 2 * k + 1 ] = c_2 + shift[ 2 * k + 1 ];
        }
        // Rprintf( "i = %d; ", i );
        // Rprintf( "idx = %d; ", ( *it_s )[ 4 * j ] );
        // Rprintf( "x = %d; ", ( *it_s )[ 4 * j + 1 ] );
        // Rprintf( "y = %d; ", ( *it_s )[ 4 * j + 2 ] );
        // Rprintf( "z = %d; ", ( *it_s )[ 4 * j + 3 ] );
        // Rprintf( "chull[0] = %d; ", chull_input[ 2 * j ] );
        // Rprintf( "chull[1] = %d;\n", chull_input[ 2 * j + 1 ] );
      }
      vector<double> hull = cHull( chull_input );
      hull_area = polygonArea( hull );
      solidity[ i ] += hull_area;
      c_spread = spread( region, ptr_seg_copy,
                         len, ptr_nidx, i );
      // if( c_spread < 1 ) {
      //   c_spread = 1;
      // }
      avg_spread[ i ] += c_spread;
      c_roundness = roundness( region, i, ptr_aidx );
      // if( c_roundness < 1 ) {
      //   c_roundness = 1;
      // }
      avg_round[ i ] += c_roundness;
      // Rprintf( "i = %d; ", i );
      // Rprintf( "slice_size = %d; ", len_slice );
      // Rprintf( "hull_area = %f; ", hull_area );
      // Rprintf( "c_spread = %f; ", c_spread );
      // Rprintf( "c_roundness =%f;\n", c_roundness );
    }
    
    solidity[ i ] /= size_region;
    avg_spread[ i ] /= p_slice.size();
    avg_round[ i ] /= p_slice.size();
    // Rprintf( "solid[%d] = %f; avg_spread[ %d ] = %f; avg_round[ %d ] = %f\n",
    //          i, solidity[ i ], i, avg_spread[ i ], i, avg_round[ i ] );
  }
}