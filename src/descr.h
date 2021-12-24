#ifndef DESCR_H
#define DESCR_H

#include <vector>

using std::vector;
// Find shape descriptors for removing regions in the last step
vector<double> descr2D( const int &len,
                        int *ptr_seg_copy, const int &label, 
                        vector<int> &region, 
                        const int *ptr_nidx, const int *ptr_aidx );
#endif