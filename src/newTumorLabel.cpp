#include "newTumorLabel.h"

// find a new tumor label for current voxel
void newTumorLabel( int &region_label, 
                    int &start, const vector<int> &tumor_labels ) {
  int i;
  for( i = start; i < tumor_labels.size(); ++ i ) {
    if( tumor_labels[ i ] == 0 ) {
      break;
    }
  }
  start = i + 1;
  region_label = - i - 4;
  return;
}