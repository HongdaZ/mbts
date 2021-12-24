#include "search.h"

#include "Rinternals.h"

void search( const int n_region, const bool init, int *label,
             const int *nidx, queue<int>& front, const int & l,
             vector<int>& region, list<int> &tumor_nbr, 
             bool &early_return ) {
  while( !front.empty() ) {
    int index = front.front();
    front.pop();
    int n_idx;
    for( int i = 0; i < 6; ++ i ) {
      n_idx = nidx[ 6 * ( index - 1 ) + i ];
      // Don't use R constants in recursive function, i.e. NA_INTEGER
      if( n_idx !=  NA_INTEGER ) { 
        
        if( label[ 2 * ( n_idx - 1 ) ] == l &&
            label[ 2 * n_idx - 1 ] != 1 ) {
          // if not initialize regions
          if( !init ) {
            // remove idx from tumor_nbr
            if( label[ 2 * n_idx - 1 ] == 2 ) {
              tumor_nbr.remove( n_idx );
              if( n_region == 0 ) {
                if( tumor_nbr.size() == 0 ) {
                  early_return = true;
                  label[ 2 * n_idx - 1 ] = 0;
                  return;
                }
              }
            }
          }
          // Rprintf( "neighbor = %d \n", n_idx );
          front.push( n_idx );
          region.push_back( n_idx );
          label[ 2 * n_idx - 1 ] = 1;
        }
      }
    }
  }
  return;
}
// 2D version of function above
void search( const int n_region, const bool init, int *label,
             const int *nidx, const int *aidx, const int &plane,
             queue<int>& front, const int & l,
             vector<int>& region, list<int> &tumor_nbr, 
             bool &early_return ) {
  const int plane_idx = aidx[ 3 * ( front.front() - 1 ) + plane ];
  int nbr_plane_idx = -1;
  while( !front.empty() ) {
    int index = front.front();
    front.pop();
    int n_idx;
    for( int i = 0; i < 6; ++ i ) {
      n_idx = nidx[ 6 * ( index - 1 ) + i ];
      // Don't use R constants in recursive function, i.e. NA_INTEGER
      if( n_idx !=  NA_INTEGER ) { 
        nbr_plane_idx = aidx[ 3 * ( n_idx - 1 ) + plane ];
        if( nbr_plane_idx == plane_idx &&
            label[ 2 * ( n_idx - 1 ) ] == l &&
            label[ 2 * n_idx - 1 ] != 1 ) {
          // if not initialize regions
          if( !init ) {
            // remove idx from tumor_nbr
            if( label[ 2 * n_idx - 1 ] == 2 ) {
              tumor_nbr.remove( n_idx );
              if( n_region == 0 ) {
                if( tumor_nbr.size() == 0 ) {
                  early_return = true;
                  label[ 2 * n_idx - 1 ] = 0;
                  return;
                }
              }
            }
          }
          // Rprintf( "neighbor = %d \n", n_idx );
          front.push( n_idx );
          region.push_back( n_idx );
          label[ 2 * n_idx - 1 ] = 1;
        }
      }
    }
  }
  return;
}