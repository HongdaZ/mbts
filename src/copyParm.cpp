#include "copyParm.h"
#include "label2col.h"

// copy estimated parameters to SEXP
void copyParm( const int &n_health,
               vector<double> &health_parm,
               vector<double> &tumor_parm,
               vector<double> &outl_parm,
               double *ptr_res_parm, const int &nrow, 
               const list<list<int>> &tumor_regions, 
               const vector<int> &outl_labels, const int &len ) {
  int j = 0;
  int curr_label;
  int cidx;
  for( int i = 0; i < n_health; ++ i ) {
    curr_label = - i - 1;
    ptr_res_parm[ nrow * j ] = curr_label;
    for( int k = 0; k < 8; ++ k ) {
      ptr_res_parm[ nrow * j + k + 1 ] = health_parm[ 8 * i + k ];
    }
    ++ j;
  }
  for( list<list<int>>::const_iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    curr_label = it->front();
    cidx = label2col( curr_label );
    ptr_res_parm[ nrow * j ] = curr_label;
    for( int k = 0; k < 8; ++ k ) {
      ptr_res_parm[ nrow * j + k + 1 ] = tumor_parm[ 8 * cidx + k ];
    }
    ++ j;
  }
  for( int i = 0; i < len; ++ i ) {
    if( outl_labels[ i ] != 0 ) {
      curr_label = i + 1;
      ptr_res_parm[ nrow * j ] = curr_label;
      for( int k = 0; k < 2; ++ k ) {
        ptr_res_parm[ nrow * j + k + 1 ] = outl_parm[ 2 * i + k ];
      }
      for( int k = 2; k < 8; ++ k ) {
        ptr_res_parm[ nrow * j + k + 1 ] = 0;
      }
      ++ j;
    }
  }
  return;
}
// health_parm: -1, -2, and -3 or -1, -2
void copyParmHealth( const int &n_health, vector<double> &health_parm,
                     double *ptr_res_parm, const int &nrow ) {
  int j = 0;
  int curr_label;
  int cidx = 0;
  for( int i = 0; i < n_health; ++ i ) {
    curr_label = - i - 1;
    ptr_res_parm[ nrow * j ] = curr_label;
    cidx = label2col( curr_label );
    for( int k = 0; k < 8; ++ k ) {
      ptr_res_parm[ nrow * j + k + 1 ] = health_parm[ 8 * cidx + k ];
    }
    ++ j;
  }
  return;
}