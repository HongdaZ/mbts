#ifndef REGION2SLICE_H
#define REGION2SLICE_H

#include <vector>
#include <list>

using std::vector;
using std::list;

vector<list<vector<int>>>
  region2slice( const list<vector<int>> &regions,
                const int &nr, const int &nc, const int &ns );

#endif