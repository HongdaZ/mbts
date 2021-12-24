#ifndef REGIONS_H
#define REGIONS_H

#include <vector>
#include <list>

using std::vector;
using std::list;

// Find connected regions
list<vector<int>> regions( int *ptr_seg, const int &len, 
                           vector<int> &region,
                           const int &label,
                           const int *ptr_nidx,
                           const int *ptr_aidx );
// 2D version of function above
vector<list<vector<int>>> regions2D( int *ptr_seg, const int &len, 
                                     vector<int> &region,
                                     const int &label,
                                     const int *ptr_nidx,
                                     const int *ptr_aidx );
// Find one region
list<vector<int>> regions( const int *ptr_seg, const int &len, 
                           const int &label,
                           const int *ptr_aidx );

#endif