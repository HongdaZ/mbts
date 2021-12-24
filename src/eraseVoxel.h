#ifndef ERASEVOXEL_H
#define ERASEVOXEL_H

#include <list>

using std::list;

void eraseVoxel( const int &idx, const int &label, 
               list<list<int>> &tumor_regions );

#endif