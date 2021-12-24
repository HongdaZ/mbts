#include "getRegion.h"
#include "clearVector.h"

// get region and region_label
void getRegion( int &region_label, vector<int> &region, 
                const vector<int> &regions_whole, 
                const vector<int> &regions_sub,
                int &start) {
  clearVector( region );
  if( start == -1 ) {
    region_label = regions_whole.front();
    for( int i = 1; i < ( regions_whole.size() - 1 ); ++ i ) {
      region.push_back( regions_whole[ i ] );
    }
    start = 0;
  } else {
    int i = start;
    region_label = regions_sub[ i ];
    while( ( ++ i ) < regions_sub.size() &&
           regions_sub[ i ] != NA_INTEGER ) {
      region.push_back( regions_sub[ i ] );
    }
    start = i + 1;
  }
  return;
}
// get region from tumor region
void getRegion( vector<int> &region, const list<int> &t_region ) {
  clearVector( region );
  for( list<int>::const_iterator it = ++ ( t_region.begin() );
       it != t_region.end(); ++ it ) {
    region.push_back( *it );
  }
}
// get healthy region
void getRegion( vector<int> &region, const int &curr_label, 
                const int *ptr_seg, const int &len ) {
  clearVector( region );
  for( int i = 1; i <= len; ++ i ) {
    if( ptr_seg[ 2 * ( i - 1 ) ] == curr_label ) {
      region.push_back( i );
    }
  }
}