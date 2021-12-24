#include <R.h>
#include <Rinternals.h>

#include "initMV.h"

// initialize yln_, yln_i and yl_
void initMV( const vector<int> &region, double *yln_, 
             double *yln_i, double *yl_, double &sum_y,
             const double *ptr_intst, const int *ptr_nidx,
             const double *ptr_nintst, const int *ptr_seg, 
             const int &curr_label ) {
  int nrow = region.size();
  int ncol = 6;
  int idx;
  int nidx;
  int nlabel;
  for( int j = 0; j < nrow; ++ j ) {
    idx = region[ j ];
    yl_[ j ] = ptr_intst[ idx - 1 ];
    sum_y += yl_[ j ];
  }
  // healthy
  if( curr_label < 0 && curr_label > - 4 ) {
    for( int i = 0; i < 6; ++ i ) {
      for( int j = 0; j < nrow; ++ j ) {
        idx = region[ j ];
        nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == curr_label ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          } 
        } else {
          yln_[ i * nrow + j ] = 0;
          yln_i[ i * nrow + j ] = 0;
        }
      }
    }
  // tumor
  } else {
    for( int i = 0; i < 6; ++ i ) {
      for( int j = 0; j < nrow; ++ j ) {
        idx = region[ j ];
        nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          nlabel = ptr_seg[ 2 * nidx - 1 ];
          if( nlabel == 3 ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
        } else {
          yln_[ i * nrow + j ] = 0;
          yln_i[ i * nrow + j ] = 0;
        }
      }
    }
  }
  return;
}
// initialize yl, yln and ylna
void initMV( double *yl, const double *yl_, double *yln,
             double *ylna, const double *yln_,
             const double *yln_i, const int &nrow, const double &mu ) {
  int ncol = 6;
  for( int j = 0; j < nrow; ++ j ) {
    yl[ j ] = yl_[ j ] - mu;
  }
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = 0; j < nrow; ++ j ) {
      if( yln_i[ i * nrow + j ] == 1 ) {
        yln[ i * nrow + j ] = yln_[ i * nrow + j ] - mu;
      } else {
        yln[ i * nrow + j ] = 0;
      }
    }
  }
  for( int i = 0; i < ( ncol / 2 ); ++ i ) {
    for( int j = 0; j < nrow; ++ j ) {
      ylna[ i * nrow + j ] = yln[ i * nrow + j ] + 
        yln[ ( ncol - 1 - i ) * nrow + j ];
    }
  }
  return;
}
// initialize yl and yln
void initMV( double *yl, const double *yl_, double *yln,
             const double *yln_,
             const double *yln_i, const int &nrow, const double &mu ) {
  int ncol = 6;
  for( int j = 0; j < nrow; ++ j ) {
    yl[ j ] = yl_[ j ] - mu;
  }
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = 0; j < nrow; ++ j ) {
      if( yln_i[ i * nrow + j ] == 1 ) {
        yln[ i * nrow + j ] = yln_[ i * nrow + j ] - mu;
      } else {
        yln[ i * nrow + j ] = 0;
      }
    }
  }
  return;
}