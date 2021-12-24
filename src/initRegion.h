#ifndef INITREGION_H
#define INITREGION_H
#include <set>
#include <map>
#include <list>
#include <vector>

using std::set;
using std::map;
using std::list;
using std::vector;

void initRegion( vector<int> &region,
                 list<list<int>> &tumor_regions, 
                 int *ptr_seg, 
                 const int *ptr_nidx, const int &len,
                 vector<int> &labels, int &n_tumor );

#endif