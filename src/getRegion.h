#ifndef GETREGION_H
#define GETREGION_H

#include <Rinternals.h>

#include <vector>
#include <list>

using std::vector;
using std::list;

// get region and region_label
void getRegion( int &region_label, vector<int> &region, 
                const vector<int> &regions_whole, 
                const vector<int> &regions_sub, 
                int &start );
// get region from tumor region
void getRegion( vector<int> &region, const list<int> &t_region );
// get healthy region
void getRegion( vector<int> &region, const int &curr_label, 
                const int *ptr_seg, const int &len );
#endif