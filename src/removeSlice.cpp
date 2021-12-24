#include "removeSlice.h"
#include "cnctRegion.h"
#include "pad2zero.h"
#include "excldRegion.h"
#include "tissueType.h"

void removeSlice( int *ptr_seg, const int &label, const int &last_slices,
                  const double &prop, const int &len,
                  vector<int> &region, 
                  const int *ptr_nidx, const int *ptr_aidx,
                  const int &nr, const int &nc, const int &ns ) {
  int slice, n, last_nonzero;
  // Find number of voxels in each slice
  vector<int> count( ns, 0 );
  for( int i = 0; i < len; ++ i ) {
    ++ count[ ptr_aidx[ 3 * i + 2 ] - 1 ];
  }
  for( int i = 0; i < ns; ++ i ) {
    if( count[ i ] != 0 ) {
      last_nonzero = i + 1;
    }
  }
  last_nonzero = last_nonzero - last_slices;
  for( int i = 0; i < len; ++ i ) {
    slice = ptr_aidx[ 3 * i + Plane::Axial ];
    if( slice > last_nonzero ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, Plane::Axial, 
                      ptr_seg, ptr_seg, 
                      label, region ) ) {
        n = region.size();
        if( n > prop * count[ slice - 1 ] ) {
          excldRegion( region, ptr_seg, 0 );
        }
      }
    }
  }
  pad2zero( ptr_seg, len );
}