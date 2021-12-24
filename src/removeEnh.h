#ifndef REMOVEENH_H
#define REMOVEENH_H

#include <vector>

using std::vector;

// Remove large regions of enh in axial planes
void removeEnh( int *ptr_seg, const int &label,
                const double &prop_brain, const double &prop_radius, 
                const double &round_rm, const int &len,
                vector<int> &region, const int *ptr_nidx, 
                const int *ptr_aidx,
                const int &nr, const int &nc, const int &ns,
                int *ptr_seg_copy );

#endif