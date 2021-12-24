#ifndef TUMORNBR_H
#define TUMORNBR_H

#include <list>

using std::list;

// find the index of the first tumor neighbor
int firstTumorNbr( const int curr_idx, const int *ptr_seg, 
              const int *ptr_nidx );
void tumorNbrLabel( list<int> &tumor_label, list<int> &tumor_nbr,
                    const int curr_idx, const int *ptr_seg,
                    const int *ptr_nidx );

#endif
