#ifndef SEARCH_H
#define SEARCH_H

#include <list> 
#include <vector>
#include <queue>

using std::list;
using std::vector;
using std::queue;

void search( const int n_region, const bool init, int *label,
             const int *nidx, queue<int>& front, const int & l,
             vector<int>& region, list<int> &tumor_nbr, 
             bool &early_return );
// 2D version of function above
void search( const int n_region, const bool init, int *label,
             const int *nidx, const int *aidx, const int &plane,
             queue<int>& front, const int & l,
             vector<int>& region, list<int> &tumor_nbr, 
             bool &early_return );

#endif
