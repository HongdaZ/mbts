#include "label2col.h"

// tranform label to the index of column
int label2col( const int label ) {
  int col;
  if( label <= - 4 ) {
    col = - label - 4;
  } else if( label < 0 && label >= - 3 ) {
    col = - label - 1;
  } else if( label > 0 ) {
    col = label - 1;
  }
  return col;
}