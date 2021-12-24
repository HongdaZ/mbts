#include "eraseOutl.h"

// erase single outlier
void eraseOutl( const int &out_label, 
                vector<int> &outl_labels,
                vector<double> &outl_parm,
                int &n_outl ) {
  int cidx = out_label - 1;
  outl_labels[ cidx ] = 0;
  for( int i = 0; i < 2; ++ i ) {
    outl_parm[ 2 * cidx + i ] = 0;
  }
  -- n_outl;
  return;
}