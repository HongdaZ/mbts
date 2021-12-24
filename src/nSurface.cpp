#include "Rinternals.h"

#include "nSurface.h"

int nSurface( int *ptr_keep, const int &len, const int *ptr_nidx ) {
  int n = 0, nidx = 0;
  bool rm;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_keep[ 2 * i ] == 1 ) {
      rm = false;
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_nidx[ 6 * i + j ];
        if( nidx != NA_INTEGER ) {
          if( ptr_keep[ 2 * ( nidx - 1 ) ] == 0 ) {
            rm = true;
            break;
          }
        } else {
          rm = true;
          break;
        }
      }
      if( rm ) {
        ++ n;
      }
    }
  }
  return n;
}