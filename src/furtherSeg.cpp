#include "furtherSeg.h"
#include "tissueType.h"

void furtherSeg( int *ptr_hgg, const int &len,
                 const int *ptr_seg, const double &prop ) {
  int n_tumor = 0, n_net = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      ++ n_tumor;
      if( ptr_seg[ 2 * i ] == Seg::SNET ) {
        ++ n_net;
      }
    }
  }
  if( n_net > ( prop * n_tumor ) ) {
    ptr_hgg[ 0 ] = -1;
  } else {
    ptr_hgg[ 0 ] = -2;
  }
}