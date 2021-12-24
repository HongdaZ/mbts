#ifndef NEWTUMORLABEL_H
#define NEWTUMORLABEL_H

#include <vector>

using std::vector;

// find an outlier label for current voxel
void newTumorLabel( int &region_label, 
                    int &start, const vector<int> &tumor_labels );

#endif