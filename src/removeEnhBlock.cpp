#include "removeEnhBlock.h"
#include "cnctRegion.h"
#include "inRegion.h"
#include "zeroVector.h"
#include "pad2zero.h"
#include "tissueType.h"
#include "pTNbr.h"
#include "radius.h"

// seg1: enh, seg2: tumor \ enh
void removeEnhBlock( int *ptr_exclude,
                     int *ptr_seg1, const int &label1,
                     int *ptr_seg2, const int &label2,
                     int *ptr_seg2_copy,
                     const double &prop1, const double &prop2, 
                     const int &len, vector<int> &region, 
                     const int *ptr_nidx, 
                     const int *ptr_aidx,
                     const int &nr, const int &nc, const int &ns ) {
  vector<int> tmp_region;
  tmp_region.reserve( len );
  int *ptr_copy_seg1 = new int[ 2 * len ]();
  int *ptr_local_seg1 = new int[ 2 * len ]();
  int *ptr_enclose = new int[ 2 * len ]();
  int idx = 0, n_in = 0, n_out = 0, n_t_nbr = 0;
  double p_t_nbr = 0, r = 0;
  const double PI = 3.141592653589793;
  
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg1[ 2 * i ] == label1 ) {
      ptr_copy_seg1[ 2 * i ] = label1;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_seg1, ptr_seg1, label1,
                    region ) ) {
      zeroVector( ptr_local_seg1, len );
      zeroVector( ptr_enclose, len );
      n_in = 0;
      n_out = 0;
      for( vector<int>::const_iterator it = region.begin();
           it != region.end(); ++ it ) {
        idx = *it;
        if( idx != 0 ) {
          ptr_local_seg1[ 2 * ( idx - 1 ) ] = label1;
          ++ n_out;
        }
      }
      inRegion( ptr_enclose, len, ptr_local_seg1, label1, 
                ptr_seg2, label2, ptr_seg2_copy,
                tmp_region, ptr_nidx,
                ptr_aidx, nr, nc, ns );
      for( int j = 0; j < len; ++ j ) {
        if( ptr_enclose[ 2 * j ] == 1 ) {
          ++ n_in;
        }
      }
      if( n_in > ( double )n_out * prop1 ) {
        for( vector<int>::const_iterator it = region.begin();
             it != region.end(); ++ it ) {
          idx = *it;
          if( idx != 0 ) {
            ptr_copy_seg1[ 2 * ( idx - 1 ) ] = 0;
          }
        }
      }
    }
  }
  pad2zero( ptr_seg1, len );
  zeroVector( ptr_enclose, len );
  inRegion2D( ptr_enclose, len, ptr_seg2, label2, ptr_copy_seg1,
              label1, ptr_seg2_copy, tmp_region, ptr_nidx, ptr_aidx,
              nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, Plane::Axial, 
                    ptr_copy_seg1, ptr_copy_seg1, 
                    label1, region ) ) {
      n_out = 0;
      n_in = 0;
      for( vector<int>::const_iterator it = region.begin();
           it != region.end(); ++ it ) {
        idx = *it;
        if( idx != 0 ) {
          ++ n_out;
          if( ptr_enclose[ 2 * ( idx - 1 ) ] == 1 ) {
            ++ n_in;
          }
        }
      }
      p_t_nbr = pTNbr2D( region, ptr_seg2, label2, ptr_nidx,
                         ptr_aidx, Plane::Axial );
      n_t_nbr = p_t_nbr * n_out;
      r = radius2D( n_out );
      if( n_in > ( 1 - prop2 ) * n_out || 
          n_t_nbr > 2 * PI * r / 8 ) {
        for( vector<int>::const_iterator it = region.begin();
             it != region.end(); ++ it ) {
          idx = *it;
          if( idx != 0 ) {
            ptr_copy_seg1[ 2 * ( idx - 1 ) ] = 0;
          }
        }
      }
    }
  }
  pad2zero( ptr_copy_seg1, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_copy_seg1[ 2 * i ] == label1 ) {
      ptr_exclude[ 2 * i ] = 1;
    }
  }
  
  delete [] ptr_copy_seg1;
  delete [] ptr_local_seg1;
  delete [] ptr_enclose;
}