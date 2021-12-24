#ifndef ADDVOXEL_H
#define ADDVOXEL_H

#include <list>

using std::list;

// add a voxel to tumor region
void addVoxel( const int &idx, const int &label, 
               list<list<int>> &tumor_regions,
               int *ptr_seg );

#endif