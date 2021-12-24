#include "Rinternals.h"

#include "enclose.h"
#include "cHull.h"
#include "inPoly.h"
#include "polygonArea.h"

// Find part of slices2 enclosed in slices1
void enclose( int *ptr_seg, const int &len, 
              const vector<list<vector<int>>> &slices1,
              const vector<list<vector<int>>> &slices2,
              vector<double> &volume,
              const int &in_sagittal, 
              const int &in_coronal, 
              const int &in_axial,
              const int &n_in ) {
  // Axial, sagittal and coronal
  vector<int> enclosed( 3 * len );
  vector<int> row{ 2, 3, 1, 3, 1, 2 };
  vector<double> shift{ -0.5, 0.5, 
                        -0.5, -0.5, 
                        0.5, -0.5, 
                        0.5, 0.5 };
  double hull_area;
  int plane1, plane2, len_slice1, len_slice2, v1, v2, idx, local_n_in;
  bool inside;
  double c_1, c_2;
  for( int i = 0; i < 3; ++ i ) {
    volume[ i ] = 0;
    const list<vector<int>> &p_slice1 = slices1[ i ];
    const list<vector<int>> &p_slice2 = slices2[ i ];
    for( list<vector<int>>::const_iterator it_s1 = p_slice1.begin();
         it_s1 != p_slice1.end(); ++ it_s1 ) {
      plane1 = ( *it_s1 )[ 1 + i ];
      
      len_slice1 = it_s1->size() / 4;
      vector<int> chull_input( 2 * len_slice1 );
      vector<double> chull_input_vol( 2 * len_slice1 * 4 );
      for( int j = 0; j < len_slice1; ++ j ) {
        
        c_1 = ( *it_s1 )[ 4 * j + row[ 2 * i ] ];
        c_2 = ( *it_s1 )[ 4 * j + row[ 2 * i + 1 ] ];
        chull_input[ 2 * j ] = c_1;
        chull_input[ 2 * j + 1 ] = c_2;
        for( int k = 0; k < 4; ++ k ) {
          chull_input_vol[ 8 * j + 2 * k ] = c_1 + shift[ 2 * k ];
          chull_input_vol[ 8 * j + 2 * k + 1 ] = c_2 + shift[ 2 * k + 1 ];
        }
      }
      vector<int> hull = cHull( chull_input );
      vector<double> hull_vol = cHull( chull_input_vol );
      hull_area = polygonArea( hull_vol );
      volume[ i ] += hull_area;
      
      for( list<vector<int>>::const_iterator it_s2 = p_slice2.begin();
           it_s2 != p_slice2.end(); ++ it_s2 ) {
        plane2 = ( *it_s2 )[ 1 + i ];
        if( plane1 == plane2 ) {
          
          len_slice2 = it_s2->size() / 4;
          vector<int> inpoly_points( 2 * len_slice2 );
          
          for( int j = 0; j < len_slice2; ++ j ) {
            inpoly_points[ 2 * j ] = ( *it_s2 )[ 4 * j + row[ 2 * i ] ];
            inpoly_points[ 2 * j + 1 ] = 
              ( *it_s2 )[ 4 * j + row[ 2 * i + 1 ] ];
          }
          
          // Close the convex hull
          v1 = hull[ 0 ];
          v2 = hull[ 1 ];
          hull.push_back( v1 );
          hull.push_back( v2 );
          vector<int> inside = inPoly( inpoly_points, hull );
          for( int j = 0; j < len_slice2; ++ j ) {
            if( inside[ j ] == 1 ) {
              idx = ( *it_s2 )[ 4 * j ];
              enclosed[ 3 * ( idx - 1 ) + i ] = 1;
            }
          }
          break;
        }
      }
    }
  }
  // Rprintf( "n_in = %d \n", n_in );
  for( int i = 0; i < len; ++ i ) {
    if( n_in == 0 ) {
      inside = true;
      if( in_sagittal == 1 ) {
        if( enclosed[ 3 * i ] == 0 ) {
          inside = false;
        }
      }
      if( in_coronal == 1 ) {
        if( enclosed[ 3 * i + 1 ] == 0 ) {
          inside = false;
        }
      }
      if( in_axial == 1 ) {
        if( enclosed[ 3 * i + 2 ] == 0 ) {
          inside = false;
        }
      }
      if( inside ) {
        ptr_seg[ 2 * i ] = 1; 
      }
    } else {
      local_n_in = 0;
      for( int j = 0; j < 3; ++ j ) {
        if( enclosed[ 3 * i + j ] == 1 ) {
          ++ local_n_in;
        }
      }
      if( local_n_in >= n_in ) {
        ptr_seg[ 2 * i ] = 1; 
      }
    }
  }
}