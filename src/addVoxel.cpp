#include "addVoxel.h"

// add a voxel to tumor region
void addVoxel( const int &idx, const int &label, 
               list<list<int>> &tumor_regions,
               int *ptr_seg ) {
  for( list<list<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    if( it->front() == label ) {
      it->push_back( idx );
      break;
    }
  }
  ptr_seg[ 2 * ( idx - 1 ) ] = label;
  return;
}
