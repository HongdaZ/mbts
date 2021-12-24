#include "Rinternals.h"

#include "areaRatio.h"

double areaRatio( const int *ptr_seg1, const int &label1,
                  const int *ptr_seg2, const int &label2,
                  const int &len, const int *ptr_ndix ) {
  double ratio = 0;
  int n_t_nbr = 0, n_crown_surface = 0, nidx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg2[ 2 * i ] == label2 ) {
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_ndix[ 6 * i + j ];
        if( nidx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( nidx - 1 ) ] != label2 &&
              ptr_seg1[ 2 * ( nidx - 1 ) ] == label1 ) {
            ++ n_t_nbr;
            break;
          } else if( ptr_seg2[ 2 * ( nidx - 1 ) ] != label2 &&
                     ptr_seg1[ 2 * ( nidx - 1 ) ] != label1 ) {
            ++ n_crown_surface;
            break;
          }
        } else {
          ++ n_crown_surface;
          break;
        }
      }
    }
  }
  ratio = ( double ) n_t_nbr / n_crown_surface;
  return ratio;
}