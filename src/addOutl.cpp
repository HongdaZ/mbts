#include "addOutl.h"
#include "label2col.h"

// add new outlier
void addOutl( int *ptr_seg, const int &idx,
              const int &new_out_label,
              const vector<double> &new_out_parm, 
              vector<int> &outl_labels,
              vector<double> &outl_parm,
              int &n_outl ) {
  ptr_seg[ 2 * ( idx - 1 ) ] = new_out_label;
  int cidx = label2col( new_out_label );
  outl_labels[ cidx ] = 1;
  for( int i = 0; i < 2; ++ i ) {
    outl_parm[ 2 * cidx + i ] = new_out_parm[ i ];
  }
  ++ n_outl;
  return;
}