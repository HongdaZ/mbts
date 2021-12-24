#ifndef SPREAD_H
#define SPREAD_H

#include <vector>

using std::vector;

// Measure the spread of the voxels in region
double spread( const vector<int> &region, int *ptr_seg_copy,
               const int &len, const int *ptr_nidx );
// 2D version of the function above
double spread( const vector<int> &region, int *ptr_seg_copy,
               const int &len, const int *ptr_nidx,
               const int &plane );
#endif