#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <list>

#include "remove.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "tissueType.h"
#include "zeroVector.h"
#include "spread.h"
#include "roundness.h"
#include "inRegion.h"
#include "removeSmall.h"

using std::max;
using std::min;
using std::list;

void remove( vector<int> &region, 
             int *ptr_seg2_copy,
             int *ptr_seg_copy,
             const int *ptr_aidx,
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis,
             int *ptr_enh, int *ptr_edema, const int &m_tumor,
             const int &m_enh, const int &m_enh_enc,
             const double &spread_factor,
             const int len, const int &nr, const int &nc, 
             const int &ns ) {
  
  int *ptr_enclose_tumor = new int[ 2 * len ]();
  int *ptr_tmp_tumor = new int[ 2 * len ]();
  vector<int> tmp_region;
  tmp_region.reserve( len );
  
  for( int i = 0; i < len; ++ i ) {
    ptr_tmp_tumor[ 2 * i ] = ptr_tumor[ 2 * i ];
  }
  int MAX = nr * nc * ns;
  int n_enc_t = 0;
  int max_size = 0;
  double max_spread_idx = 0, min_spread_idx = 0;
  int n_enh = 0;
  int idx = 0;
  double spread_idx = 0;
  int size = 0;
  double s_factor;
  list<int> region_size;
  list<double> region_spread;
  
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      spread_idx = spread( region, ptr_seg_copy, 
                           len, ptr_nidx );
      region_size.push_back( region.size() );
      if( max_size < region.size() ) {
        max_size = region.size();
      }
      region_spread.push_back( spread_idx );
      // min_spread is also updated later
      if( max_spread_idx < spread_idx ) {
        max_spread_idx = spread_idx;
      }
    }
  }
  pad2zero( ptr_tumor, len );
  size = min( m_tumor, max_size );
  min_spread_idx = max_spread_idx;
  
  list<int>::const_iterator it_size = region_size.begin();
  list<double>::const_iterator it_spread = region_spread.begin();
  
  for( ; it_size != region_size.end(); ++ it_size, ++ it_spread ) {
    if( *it_size > 5000 ) {
      spread_idx = *it_spread;
      if( spread_idx < min_spread_idx ) {
        min_spread_idx = spread_idx;
      }
    }
  }
  
  s_factor = max( min_spread_idx, spread_factor );
  // Rprintf( "s_factor = %f\n", s_factor );
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      spread_idx = spread( region, ptr_seg_copy,
                           len, ptr_nidx );
      // Rprintf( "spread_idx = %f\n", spread_idx );
      if( region.size() < size &&
          spread_idx > s_factor ) {
        excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                     ptr_necrosis, ptr_enh, ptr_edema, MAX );
      }
    }
  }
  pad2zero( ptr_tumor, len );
  // Tumor enclosed by enh (including enh)
  inRegion( ptr_enclose_tumor, len, ptr_enh, Tumor::ET,
            ptr_tmp_tumor, 1, ptr_seg2_copy,
            tmp_region, ptr_nidx, ptr_aidx,
            nr, nc, ns );
  // Rprintf( "s_factor = %f\n", s_factor );
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      if( region.size() < size ) {
        // Find enh inside the tumor region
        n_enh = 0;
        for( int j = 0; j < region.size(); ++ j ) {
          idx = region[ j ];
          if( idx != 0 ) {
            if( ptr_enh[ 2 * ( idx - 1 ) ] == Tumor::ET ) {
              ++ n_enh;
            }
          }
        }
        if( n_enh < m_enh ) {
          // If not including enh, exclude small sized regions
          if( region.size() < size ) {
            excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                         ptr_necrosis, ptr_enh, ptr_edema, MAX );
          }
        } else {
          // Remove regions with tumor enclosed enh < m_enh_enc
          n_enc_t = 0;
          for( int j = 0; j < region.size(); ++ j ) {
            idx = region[ j ];
            if( idx != 0 ) {
              if( ptr_enclose_tumor[ 2 * ( idx - 1 ) ] == 1 ) {
                ++ n_enc_t;
              }
            }
          }
          if( n_enc_t < m_enh_enc ) {
            excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                         ptr_necrosis, ptr_enh, ptr_edema, MAX );
          }
        }
      }
    }
  }
  pad2zero( ptr_tumor, len );
  
  delete [] ptr_enclose_tumor;
  delete [] ptr_tmp_tumor;
}
// Remove regions with size < m_tumor &&
//                     ( spread_idx > 2 ||
//                       ( max enh size < m_enh &&
//                         enclosed tumor < m_enh_enc &&
//                         prop enclose tumor < 0.75 ) )
void remove( vector<int> &region, 
             int *ptr_seg2_copy,
             int *ptr_seg_copy,
             const int *ptr_aidx,
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis,
             int *ptr_enh, int *ptr_edema, const int &m_tumor,
             const int &m_enh, const int &m_enh_enc, const int len,
             const int &nr, const int &nc, const int &ns,
             const double &rm_spread, const double &rm_round ) {
  
  int *ptr_enclose_tumor = new int[ 2 * len ]();
  int *ptr_tmp_tumor = new int[ 2 * len ]();
  vector<int> tmp_region;
  tmp_region.reserve( len );
  int max_size = 0, n_enc_t = 0, idx = 0;
  int MAX = nr * nc * ns;
  bool remove = true;
  double spread_idx = 0, round_idx = 0;
  
  for( int i = 0; i < len; ++ i ) {
    ptr_tmp_tumor[ 2 * i ] = ptr_tumor[ 2 * i ];
  }
  
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tmp_tumor, ptr_tmp_tumor,
                    1, region ) ) {
      if( max_size < region.size() ) {
        max_size = region.size();
      }
    }
  }
  pad2zero( ptr_tmp_tumor, len );
  inRegion( ptr_enclose_tumor, len, ptr_enh, Tumor::ET,
            ptr_tmp_tumor, 1, ptr_seg2_copy,
            tmp_region, ptr_nidx, ptr_aidx,
            nr, nc, ns );
  max_size = min( m_tumor, max_size );
  // remove
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      if( region.size() < max_size ) {
        remove = true;
        if( region.size() > 200 ) {
          spread_idx = spread( region, ptr_seg_copy,
                               len, ptr_nidx );
          round_idx = roundness( region, ptr_aidx );
          if( spread_idx < rm_spread  && 
              round_idx < rm_round ) {
            // Find enh inside the tumor region
            for( int j = 0; j < region.size(); ++ j ) {
              idx = region[ j ];
              if( idx != 0 ) {
                cnctRegion( idx, ptr_nidx, ptr_enh, ptr_enh, Tumor::ET,
                            tmp_region );
                if( tmp_region.size() > m_enh ) {
                  remove = false;
                  break;
                }
              }
            }
            pad2zero( ptr_enh, len );
            if( remove ) {
              for( int j = 0; j < region.size(); ++ j ) {
                idx = region[ j ];
                if( idx != 0 ) {
                  cnctRegion( idx, ptr_nidx,
                              ptr_enclose_tumor, ptr_enclose_tumor, 1,
                              tmp_region );
                  if( tmp_region.size() > m_enh_enc ) {
                    remove = false;
                    break;
                  } else if( ( double ) tmp_region.size() /
                    region.size() > 0.80 ){
                    remove = false;
                    break;
                  }
                }
              }
              pad2zero( ptr_enclose_tumor, len );
            }
          } 
        }
        if( remove ) {
          excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                       ptr_necrosis, ptr_enh, ptr_edema, MAX );
        }
      }
    }
  }
  pad2zero( ptr_tumor, len );
  
  delete [] ptr_enclose_tumor;
  delete [] ptr_tmp_tumor;
  
}