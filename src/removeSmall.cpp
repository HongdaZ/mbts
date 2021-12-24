#include "removeSmall.h"
#include "excldRegion.h"
#include "cnctRegion.h"
#include "pad2zero.h"

void removeSmall( vector<int> &region,
                  const int *ptr_nidx,
                  int *ptr_seg, int *ptr_tumor,
                  int *ptr_hemorrhage, int *ptr_necrosis,
                  int *ptr_enh, int *ptr_edema, const int &size, 
                  const int &len ) {
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                   ptr_necrosis, ptr_enh, ptr_edema, size );
    }
  }
  pad2zero( ptr_tumor, len );
}