#ifndef FINDREGION_H
#define FINDREGION_H

#include <vector>
#include <list>
#include <queue>

using std::list;
using std::vector;
using std::queue;

void findRegion( vector<int> &region,
                 const int &n_region, int *ptr_label, const int *ptr_nidx, 
                 const bool &init, list<int> &tumor_nbr, bool &early_return, 
                 const int &start );
// 2D version of function above
void findRegion( vector<int> &region,
                 const int &n_region, int *ptr_label, const int *ptr_nidx, 
                 const int *ptr_aidx, const int &plane,
                 const bool &init, list<int> &tumor_nbr, bool &early_return, 
                 const int &start );

#endif