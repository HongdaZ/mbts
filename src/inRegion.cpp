#include <list>

#include "inRegion.h"
#include "regions.h"
#include "region2slice.h"
#include "enclose.h"
#include "getRange.h"
#include "zeroVector.h"
#include "smallRegion.h"

using std::list;
using std::vector;

// Find the points of ptr_seg2 enclosed in ptr_seg1
vector<double> inRegion( int *ptr_enclose, const int &len,
                         int *ptr_seg1, const int &label1, 
                         int *ptr_seg2, const int &label2,
                         int *ptr_seg2_copy,
                         vector<int> &region, 
                         const int *ptr_nidx, const int *ptr_aidx,
                         const int &nr, const int &nc, const int &ns,
                         const int &in_sagittal, 
                         const int &in_coronal, 
                         const int &in_axial,
                         const int &n_in ) { 
  vector<double> volume( 3, 0 );
  // Find connected regions in ptr_seg1
  list<vector<int>> regions1 = regions( ptr_seg1, len,
                                        region, label1, 
                                        ptr_nidx,
                                        ptr_aidx );
  vector<list<vector<int>>> slices1 =
    region2slice( regions1, nr, nc, ns );
  
  zeroVector( ptr_seg2_copy, len );
  vector<int> range = getRange( ptr_seg1, label1, ptr_aidx, len,
                                nr, nc, ns );
  smallRegion( range, ptr_seg2, label2, ptr_aidx, ptr_seg2_copy,
               len );
  
  list<vector<int>> regions2 = regions( ptr_seg2_copy, len,
                                        label2, ptr_aidx ); 
  vector<list<vector<int>>> slices2 = 
    region2slice( regions2, nr, nc, ns );
  enclose( ptr_enclose, len, slices1, slices2, volume,
           in_sagittal, in_coronal, in_axial, n_in );
  return volume;
}
// 2D version of inregion
vector<double> inRegion2D( int *ptr_enclose, const int &len,
                           int *ptr_seg1, const int &label1, 
                           int *ptr_seg2, const int &label2,
                           int *ptr_seg2_copy,
                           vector<int> &region, 
                           const int *ptr_nidx, const int *ptr_aidx,
                           const int &nr, const int &nc, const int &ns,
                           const int &in_sagittal, 
                           const int &in_coronal, 
                           const int &in_axial,
                           const int &n_in ) {
  vector<double> volume( 3, 0 );
  // Find connected regions in ptr_seg1
  vector<list<vector<int>>> slices1 = regions2D( ptr_seg1, len,
                                        region, label1, 
                                        ptr_nidx,
                                        ptr_aidx );
  zeroVector( ptr_seg2_copy, len );
  vector<int> range = getRange( ptr_seg1, label1, ptr_aidx, len,
                                nr, nc, ns );
  smallRegion( range, ptr_seg2, label2, ptr_aidx, ptr_seg2_copy,
               len );
  
  list<vector<int>> regions2 = regions( ptr_seg2_copy, len,
                                        label2, ptr_aidx ); 
  vector<list<vector<int>>> slices2 = 
    region2slice( regions2, nr, nc, ns );
  enclose( ptr_enclose, len, slices1, slices2, volume,
           in_sagittal, in_coronal, in_axial, n_in );
  return volume;
}