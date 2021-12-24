#ifndef LABELREGION_H
#define LABELREGION_H

#include <vector>

using std::vector;

void labelRegion( const vector<int> &region,
                  int *ptr_seg );
void recoverLabel( const vector<int> &region, 
                   int *ptr_seg );


#endif