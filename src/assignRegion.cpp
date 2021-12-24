#include <Rinternals.h>

#include "assignRegion.h"


// store region in whole and sub-regions
void assignRegion( vector<int> &regions_whole,  
                   vector<int> &regions_sub,
                   const vector<int> &region,
                   const int &label ) {
  regions_sub.push_back( label );
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    regions_whole.push_back( *it );
    regions_sub.push_back( *it );
  }
  regions_sub.push_back( NA_INTEGER );
  return;
}

void assignRegion( vector<int> &regions_whole,  
                   vector<int> &regions_sub,
                   const list<list<int>> &tumor_regions,
                   const int &label,
                   int &r_size ) {
  list<list<int>>::const_iterator it = tumor_regions.begin();
  while( it != tumor_regions.end() && 
         it->front() != label ) {
    ++ it;
  }
  list<int>::const_iterator it_region = ++ ( it->begin() );
  regions_sub.push_back( label );
  r_size = 0;
  for( ; it_region != it->end(); ++ it_region ) {
    regions_whole.push_back( *it_region );
    regions_sub.push_back( *it_region );
    ++ r_size;
  }
  regions_sub.push_back( NA_INTEGER );
  return;
}