#ifndef ADDINSIDE_H
#define ADDINSIDE_H

#include <vector>

using std::vector;

// Add voxels of ptr_seg2 inside ptr_seg1
void addInside( int *ptr_inside, const int &len,
                int *ptr_seg1, const int &label1, 
                int *ptr_seg2, const int &label2,
                int *ptr_seg2_copy,
                vector<int> &region, 
                const int *ptr_nidx, const int *ptr_aidx,
                const int &nr, const int &nc, const int &ns );

#endif