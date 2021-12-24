#ifndef ROUNDNESS_H
#define ROUNDNESS_H

#include <vector>

using std::vector;

// Measure the roundness of the voxels in region
double roundness( const vector<int> &region, const int *ptr_aidx );
// 2D version of the function above
double roundness( const vector<int> &region, const int &plane,
                  const int *ptr_aidx );
#endif