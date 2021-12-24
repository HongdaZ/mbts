#include "excldRegion.h"

#include "Rinternals.h"

int excldRegion( const vector<int> &region, const int *ptr_nidx,
                  int *ptr_seg2,
                  const int *ptr_seg1, const int &label1, 
                  const bool &has ) {
  int n_idx, index;
  bool remove = has;
  
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 ) {
      for( int i = 0; i < 6; ++ i ) {
        n_idx = ptr_nidx[ 6 * ( index - 1 ) + i ];
        if( n_idx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( n_idx - 1 ) ] == 0 &&
              ptr_seg1[ 2 * ( n_idx - 1 ) ] == label1 ) {
            if( has ) {
              return n_idx;
            } else {
              remove = true;
              break;
            }
          }
        }
      }
    }
  }
  if( remove ) {
    for( vector<int>::const_iterator it = region.begin(); 
         it != region.end(); ++ it ) {
      index = *it;
      if( index != 0 ) {
        ptr_seg2[ 2 * ( index - 1 ) ] = 0;
      }
    }
    return 0;
  } else {
    return -1;
  }
  
}
int excldRegion( const vector<int> &region, const int *ptr_nidx,
                 const int *ptr_seg1, const int &label1, int *ptr_seg2 ) {
  int n_idx, index;
  bool remove = false;
  
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 ) {
      for( int i = 0; i < 6; ++ i ) {
        n_idx = ptr_nidx[ 6 * ( index - 1 ) + i ];
        if( n_idx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( n_idx - 1 ) ] == 0 &&
              ptr_seg1[ 2 * ( n_idx - 1 ) ] == 0 ) {
            remove = true;
            break;
          }
        }
      }
    }
  }
  if( remove ) {
    for( vector<int>::const_iterator it = region.begin(); 
         it != region.end(); ++ it ) {
      index = *it;
      if( index != 0 ) {
        ptr_seg2[ 2 * ( index - 1 ) ] = 0;
      }
    }
    return 0;
  } else {
    return -1;
  }
  
}
void excldRegion( const vector<int> &region,
                  int *ptr_seg1,
                  int *ptr_seg2, 
                  const int &label, const int &size, const double &prop,
                  int *ptr_hemorrhage,
                  int *ptr_necrosis,
                  int *ptr_enh,
                  int *ptr_edema ) {
  int index, n = 0;
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 && ptr_seg2[ 2 * ( index - 1 ) ] == label ) {
      ++ n;
    }
  }
  if( n > size && ( n / ( double )region.size() ) < prop ) {
    return;
  }
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 ) {
      ptr_seg1[ 2 * ( index - 1 ) ] = 0;
      ptr_seg2[ 2 * ( index - 1 ) ] = 0;
      ptr_hemorrhage[ 2 * ( index - 1 ) ] = 0;
      ptr_necrosis[ 2 * ( index - 1 ) ] = 0;
      ptr_enh[ 2 * ( index - 1 ) ] = 0;
      ptr_edema[ 2 * ( index - 1 ) ] = 0;
    }
  }
}

void excldRegion( const vector<int> &region,
                  int *ptr_seg, int *ptr_tumor,
                  int *ptr_hemorrhage, int *ptr_necrosis,
                  int *ptr_enh, int *ptr_edema, const int &size ) {
  if( region.size() > size ) {
    return;
  } 
  int index;
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 ) {
      ptr_seg[ 2 * ( index - 1 ) ] = 0;
      ptr_tumor[ 2 * ( index - 1 ) ] = 0;
      ptr_hemorrhage[ 2 * ( index - 1 ) ] = 0;
      ptr_necrosis[ 2 * ( index - 1 ) ] = 0;
      ptr_enh[ 2 * ( index - 1 ) ] = 0;
      ptr_edema[ 2 * ( index - 1 ) ] = 0;
    }
  }
}
void excldRegion( const vector<int> &region,
                  int *ptr_enclose_ncr, const int &size ) {
  if( region.size() < size ) {
    return;
  }
  int index;
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    index = *it;
    if( index != 0 ) {
      ptr_enclose_ncr[ 2 * ( index - 1 ) ] = 0;
    }
  }
}