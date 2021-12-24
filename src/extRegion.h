#ifndef EXTREGION_H
#define EXTREGION_H

#include <vector>

using std::vector;

void extRegion( const vector<int> &region, int *ptr_seg,
                const int &label, const double &ratio,
                bool remove = false );

#endif