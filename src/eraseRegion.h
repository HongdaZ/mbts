#ifndef ERASEREGION_H
#define ERASEREGION_H

#include <list>
#include <map>
#include <vector>

using std::list;
using std::map;
using std::vector;

// erase single tumor region
void eraseRegion( const int &tumor_label, 
                  vector<int> &tumor_labels,
                  vector<double> &tumor_parm,
                  list<list<int>> &tumor_regions,
                  int &n_tumor );

#endif