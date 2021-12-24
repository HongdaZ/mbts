#ifndef ASSIGNREGION_H
#define ASSIGNREGION_H

#include <list>
#include <vector>

using std::vector;
using std::list;

// store region in whole and sub-regions
void assignRegion( vector<int> &regions_whole,  
                   vector<int> &regions_sub,
                   const vector<int> &region,
                   const int &label );
void assignRegion( vector<int> &regions_whole,  
                   vector<int> &regions_sub,
                   const list<list<int>> &tumor_regions,
                   const int &label, 
                   int &r_size );

#endif