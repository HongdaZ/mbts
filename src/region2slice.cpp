#include "region2slice.h"

// Convert regions to slices
vector<list<vector<int>>>
  region2slice( const list<vector<int>> &regions,
                const int &nr, const int &nc, const int &ns ) {
    vector<list<vector<int>>> res( 3 );
    vector<int> dim{ nr, nc, ns };
    int np_r;
    int np_s;
    // Axial, sagittal and coronal
    for( int i = 0; i < 3; ++ i ) {
      list<vector<int>> plane;
      for( list<vector<int>>::const_iterator it_r = regions.begin();
           it_r != regions.end(); ++ it_r ) {
        np_r = it_r->size() / 4;
        for( int j = 0; j < dim[ i ]; ++ j ) {
          list<int> slice;
          for( int k = 0; k < np_r; ++ k ) {
            if( ( *it_r )[ 4 * k + 1 + i ] - 1 == j ) {
              slice.push_back( k );
            }
          }
          if( slice.size() > 0 ) {
            np_s = 0;
            vector<int> slice_v( slice.size() * 4 );
            for( list<int>::const_iterator it_s = slice.begin();
                 it_s != slice.end(); ++ it_s, ++ np_s ) {
              for( int k = 0; k < 4; ++ k ) {
                slice_v[ 4 * np_s + k ] = 
                  ( *it_r )[ 4 * *it_s + k ];
              }
            }
            plane.push_back( slice_v );
          }
        }
      }
      res[ i ] = plane;
    }
    return res;
  }