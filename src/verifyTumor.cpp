#include "verifyTumor.h"

bool verifyTumor( map<int, list<int>> &tumor_regions, const int *ptr_seg,
                  const int len ) {
  bool res = true;
  int label;
  int n = 0;
  for( map<int, list<int>>::iterator it = tumor_regions.begin(); 
       it != tumor_regions.end(); ++ it ) {
    label = it->first;
    for( list<int>::iterator it_list = it->second.begin(); 
         it_list != it->second.end(); ++ it_list ) {
      ++ n;
      if( ptr_seg[ 2 * ( *it_list - 1 ) ] != label ) {
        res = false;
        return res;
      }
    }
  }
  int m = 0;
  for( int i = 1; i <= len; ++ i ) {
    if( ptr_seg[ 2 * ( i - 1 ) ] < - 3 ) {
      ++ m;
    }
  }
  res = res && ( m == n );
  return res;
}