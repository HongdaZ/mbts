#ifndef REMOVESLICE_H
#define REMOVESLICE_H

#include <vector>

using std::vector;

// Remove large regions of enh in axial planes
void removeSlice( int *ptr_seg, const int &label, const int &last_slices,
                  const double &prop, const int &len,
                  vector<int> &region, 
                  const int *ptr_nidx, const int *ptr_aidx,
                  const int &nr, const int &nc, const int &ns );

#endif