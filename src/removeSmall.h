#ifndef REMOVESMALL_H
#define REMOVESMALL_H

#include <vector>

using std::vector;

// Remove small regions with size < size
void removeSmall( vector<int> &region,
                  const int *ptr_nidx,
                  int *ptr_seg, int *ptr_tumor,
                  int *ptr_hemorrhage, int *ptr_necrosis,
                  int *ptr_enh, int *ptr_edema, const int &size, 
                  const int &len );

#endif