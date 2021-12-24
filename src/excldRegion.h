#ifndef EXCLDREGION_H
#define EXCLDREGION_H

#include <vector>

using std::vector;

// Remove the region from ptr_seg2 if there is no adjacent voxel with 
// ptr_seg1 == label1
int excldRegion( const vector<int> &region, const int *ptr_nidx,
                 int *ptr_seg2,
                 const int *ptr_seg1, 
                 const int &label1, const bool &has = true );
// Remove the region from ptr_seg2 if there are neighboring voxels with 
// ptr_seg1 == 0 && ptr_seg2 == 0
int excldRegion( const vector<int> &region, const int *ptr_nidx,
                 const int *ptr_seg1, const int &label1, int *ptr_seg2 );
// Remove the region from ptr_seg1 if the number of voxels with 
// ptr_seg2 == label is less than size
void excldRegion( const vector<int> &region,
                  int *ptr_seg1,
                  int *ptr_seg2, 
                  const int &label, const int &size, const double &prop,
                  int *ptr_hemorrhage,
                  int *ptr_necrosis,
                  int *ptr_enh,
                  int *ptr_edema ) ;
// Remove tumor regions with size < size.
void excldRegion( const vector<int> &region,
                  int *ptr_seg, int *ptr_tumor,
                  int *ptr_hemorrhage, int *ptr_necrosis,
                  int *ptr_enh, int *ptr_edema, const int &size );
// Remove newly added necrosis if the region size > size
void excldRegion( const vector<int> &region,
                  int *ptr_enclose_ncr, const int &size );

#endif