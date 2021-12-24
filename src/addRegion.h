#ifndef ADDREGION_H
#define ADDREGION_H

#include <Rinternals.h>

#include <list>
#include <vector>

using std::list;
using std::vector;

// add new tumor regions
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const int &n_row, 
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                list<list<int>> &tumor_regions,
                const vector<int> &regions_sub,
                int &n_tumor );
// add a new single voxel tumor region
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const int &idx,
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                list<list<int>> &tumor_regions,
                int &n_tumor );
#endif