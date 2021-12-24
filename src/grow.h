#ifndef GROW_H
#define GROW_H

#include <vector>

using std::vector;
// Let the apple grow
void grow( int *ptr_seg_copy,
           const int &last, int &n_tumor, const int &len, 
           vector<int> &region,
           const int *ptr_nidx, const int *ptr_aidx,
           int *ptr_whole, int *ptr_res,
           int *ptr_one, int *ptr_keep, int *ptr_remain, 
           const double &s_trim, const double &r_trim ); 

#endif