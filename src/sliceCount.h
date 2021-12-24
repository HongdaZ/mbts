#ifndef SLICECOUNT_H
#define SLICECOUNT_H

#include <vector>

using std::vector;

void sliceCount( vector<int> &counts, 
                 int *ptr_seg, const int &len, const int &label,
                 const int *ptr_nidx, const int *ptr_aidx,
                 const int &nr, const int &nc, const int &ns );
void sliceCount( vector<int> &counts, const vector<int> &region,
                 const int *ptr_nidx, const int *ptr_aidx,
                 const int &nr, const int &nc, const int &ns );

#endif