#ifndef NBRLABEL_H
#define NBRLABEL_H

#include <list>

using std::list;

void nbrLabel( list<int> &nbr_label, list<int> &tumor_nbr, 
               const int curr_idx, const int *ptr_seg, 
               const int *ptr_nidx );
                 
#endif
