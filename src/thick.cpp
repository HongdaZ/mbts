#include "Rinternals.h"
#include <list>

#include "thick.h"
#include "clearVector.h"
#include "pad2zero.h"

using std::list;

int thick( const int *ptr_seg1, const int &label1, 
           int *ptr_seg2, const int &label2,
           const int &len, const int *ptr_nidx ) {
  int thickness = 0;
  list<int> outer_old, outer_new;
  int idx = 0, nidx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg2[ 2 * i ] == label2 ) {
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_nidx[ 6 * i + j ];
        if( nidx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( nidx - 1 ) ] != label2 &&
              ptr_seg1[ 2 * ( nidx - 1 ) ] == label1 ) {
            outer_old.push_back( i + 1 );
            ptr_seg2[ 2 * i + 1 ] = 1; 
            break;
          }
        }
      }
    }
  }
  while( outer_old.size() > 0 ) {
    ++ thickness;
    for( list<int>::const_iterator it = outer_old.begin();
         it != outer_old.end(); ++ it ) {
      idx = *it;
      for( int i = 0; i < 6; ++ i ) {
        nidx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( nidx - 1 ) ] == label2 &&
              ptr_seg2[ 2 * nidx - 1 ] != 1 ) {
            outer_new.push_back( nidx );
            ptr_seg2[ 2 * nidx - 1 ] = 1;
          }
        }
      }
    }
    outer_old = outer_new;
    clearVector( outer_new );
  }
  pad2zero( ptr_seg2, len );
  return thickness;
}