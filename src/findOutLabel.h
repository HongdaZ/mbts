#ifndef FINDOUTLABEL_H
#define FINDOUTLABEL_H

#include <vector>

using std::vector;

// find an outlier label for current voxel
void findOutLabel( int &out_label, const int &idx, const int *ptr_seg, 
                   const vector<int> &outl_labels );

#endif