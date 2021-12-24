#ifndef ONREGION_H
#define ONREGION_H

#include <vector>

using std::vector;

void onRegion( int *ptr_on, const int &len, const double &prop,
               int *ptr_seg1, const int &label1, 
               int *ptr_seg2, const int &label2,
               int *ptr_seg2_copy,
               vector<int> &region, const double &spread_factor,
               const int *ptr_nidx, const int *ptr_aidx,
               const int &nr, const int &nc, const int &ns,
               int *ptr_seg_copy, const double &prop_hull = 0.3,
               const double &nt_prop = 0.3 );

#endif