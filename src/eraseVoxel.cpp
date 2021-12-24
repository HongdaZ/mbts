#include "eraseVoxel.h"

// remove a  voxel to tumor region
void eraseVoxel( const int &idx, const int &label, 
               list<list<int>> &tumor_regions ) {
  for( list<list<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    if( it->front() == label ) {
      list<int>::iterator it_region = it->begin();
      while( ( ++ it_region ) != it->end() ) {
        if( *it_region == idx ) {
          it->erase( it_region );
          break;
        }
      }
      break;
    }
  }
  return;
}
