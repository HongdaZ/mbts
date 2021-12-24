#include "addRegion.h"
#include "label2col.h"

// add new tumor regions
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const int &n_row, 
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                list<list<int>> &tumor_regions,
                const vector<int> &regions_sub,
                int &n_tumor ) {
  int ncol = region_parm.size() / n_row;
  int new_label;
  int cidx;
  int start = 0;
  for( int i = 0; i < ncol; ++ i ) {
    list<int> new_region;
    new_label = region_parm[ n_row * i ];
    new_region.push_back( new_label );
    cidx = label2col( new_label );
    tumor_labels[ cidx ] = 1;
    for( int j = 0; j < 8; ++ j ) {
      tumor_parm[ 8 * cidx + j ] = region_parm[ n_row * i + j + 1 ];
    }
    // start at label, start + 1 at idx
    for( int j = start + 1; j < regions_sub.size(); ++ j ) {
      if( regions_sub[ j ] == NA_INTEGER ) {
        tumor_regions.push_back( new_region );
        // j at NA_INTEGER, j + 1 at label
        start = j + 1;
        break;
      } else {
        new_region.push_back( regions_sub[ j ] );
        ptr_seg[ 2 * ( regions_sub[ j ] - 1 ) ] = new_label;
      }
    }
    ++ n_tumor;
  }
  return;
}
// add a new single voxel tumor region
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const int &idx,
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                list<list<int>> &tumor_regions,
                int &n_tumor ) {
  int new_label = region_parm[ 0 ];
  int cidx = label2col( new_label );
  tumor_labels[ cidx ] = 1;
  for( int j = 0; j < 8; ++ j ) {
    tumor_parm[ 8 * cidx + j ] = region_parm[ j + 1 ];
  }
  ++ n_tumor;
  list<int> new_region;
  new_region.push_back( new_label );
  new_region.push_back( idx );
  tumor_regions.push_back( new_region );
  ptr_seg[ 2 * ( idx - 1 ) ] = new_label;
  return;
}