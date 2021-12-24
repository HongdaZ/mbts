#ifndef GETRANGE_H
#define GETRANGE_H

#include <vector>

using std::vector;

vector<int> getRange( const int *ptr_seg, const int &label,
                      const int *ptr_aidx, 
                      const int &len,
                      const int &nr, const int &nc, const int &ns );

#endif