#include "pad2zero.h"

// Recover the padding to zero's
void pad2zero( int *ptr_seg, const int &len ) {
  for( int i = 0; i < len; ++ i ) {
    ptr_seg[ 2 * i + 1 ] = 0;
  }
}
void pad2zero( int *ptr_seg, const vector<int> region ) {
  int idx = 0;
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ptr_seg[ 2 * idx - 1 ] = 0;
    }
  }
}