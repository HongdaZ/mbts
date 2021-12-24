#ifndef REMOVEENHBLOCK_H
#define REMOVEENHBLOCK_H

#include <vector>

using std::vector;

// Remove large blocks of ptr_seg1
void removeEnhBlock( int *ptr_exclude,
                     int *ptr_seg1, const int &label1,
                     int *ptr_seg2, const int &label2,
                     int *ptr_seg2_copy,
                     const double &prop1, const double &prop2, 
                     const int &len, vector<int> &region, 
                     const int *ptr_nidx, 
                     const int *ptr_aidx,
                     const int &nr, const int &nc, const int &ns );

#endif