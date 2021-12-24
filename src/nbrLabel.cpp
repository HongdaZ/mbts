#include <R.h>
#include <Rinternals.h>

#include <list>
#include "nbrLabel.h"

using std::list;
// find the non-nan labels of the neighbors and indeces of 
// tumor neighbors
void nbrLabel( list<int> &nbr_label, list<int> &tumor_nbr,
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
      }
      nbr_label.push_back( label );
    }
  }
  nbr_label.sort();
  nbr_label.unique();
  return;
}
