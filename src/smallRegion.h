#ifndef SMALLREGION_H
#define SMALLREGION_H

#include <vector>

using std::vector;

void smallRegion( const vector<int> &range, 
                  const int *ptr_seg2, const int &label2, 
                  const int *ptr_aidx,
                  int *ptr_seg2_copy, const int &len );

#endif