#ifndef CNCTREGION_H
#define CNCTREGION_H

#include <vector>

using std::vector;

bool cnctRegion( const int &start, const int *ptr_nidx,
                 int *ptr_seg1, int *ptr_seg2, 
                 const int &label, vector<int> &region );
// 2D version of function above
bool cnctRegion( const int &start, const int *ptr_nidx, 
                 const int *ptr_aidx, const int &plane,
                 int *ptr_seg1, int *ptr_seg2, 
                 const int &label, vector<int> &region );
#endif