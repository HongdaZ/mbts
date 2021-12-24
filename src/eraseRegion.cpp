#include "eraseRegion.h"
#include "label2col.h"

// erase single tumor region
void eraseRegion( const int &tumor_label, 
                  vector<int> &tumor_labels,
                  vector<double> &tumor_parm,
                  list<list<int>> &tumor_regions,
                  int &n_tumor ) {
  int cidx = label2col( tumor_label );
  tumor_labels[ cidx ] = 0;
  for( int i = 0; i < 8; ++ i ) {
    tumor_parm[ 8 * cidx + i ] = 0;
  }
  for( list<list<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    if( it->front() == tumor_label ) {
      tumor_regions.erase( it );
      break;
    }
  }
  -- n_tumor;
  return;
}