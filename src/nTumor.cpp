#include "nTumor.h"

int nTumor( int *ptr_whole, const int &len ) {
  int n = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_whole[ 2 * i ] == 1 ) {
      ++ n;
    }
  }
  return n;
}