#include <R.h>
#include <Rinternals.h>

#include "tumorNbr.h"
// find the index of the first tumor neighbor
int firstTumorNbr( const int curr_idx, const int *ptr_seg, 
              const int *ptr_nidx ) {
  int nbr_idx = 0;
  int label = 0;
  int first_tumor_nbr = 0;
  for( int i = 0; i < 6; ++ i ) {
    nbr_idx = ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nbr_idx != NA_INTEGER ) {
      label = ptr_seg[ 2 * ( nbr_idx - 1 ) ];
      if( label <= - 4 ) {
        first_tumor_nbr = nbr_idx;
        break;
      }
    }
  }
  return first_tumor_nbr;
}
void tumorNbrLabel( list<int> &tumor_label, list<int> &tumor_nbr,
               const int curr_idx, const int *ptr_seg,
               const int *ptr_nidx ) {
  int nbr_idx = 0;
  int label;
  for( int i = 0; i < 6; ++ i ) {
    nbr_idx = ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nbr_idx != NA_INTEGER ) {
      label = ptr_seg[ 2 * ( nbr_idx - 1 ) ];
      if( label <= - 4 ) {
        tumor_nbr.push_back( nbr_idx );
        tumor_label.push_back( label );
      }
      
    }
  }
  tumor_label.sort();
  tumor_label.unique();
  return;
}
