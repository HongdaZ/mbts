#ifndef INREGION_H
#define INREGION_H

#include <vector>

using std::vector;

vector<double> inRegion( int *ptr_enclose, const int &len,
                         int *ptr_seg1, const int &label1, 
                         int *ptr_seg2, const int &label2,
                         int *ptr_seg2_copy,
                         vector<int> &region, 
                         const int *ptr_nidx, const int *ptr_aidx,
                         const int &nr, const int &nc, const int &ns,
                         const int &in_sagittal = 0, 
                         const int &in_coronal = 0, 
                         const int &in_axial = 0,
                         const int &n_in = 3 );
vector<double> inRegion2D( int *ptr_enclose, const int &len,
                           int *ptr_seg1, const int &label1, 
                           int *ptr_seg2, const int &label2,
                           int *ptr_seg2_copy,
                           vector<int> &region, 
                           const int *ptr_nidx, const int *ptr_aidx,
                           const int &nr, const int &nc, const int &ns,
                           const int &in_sagittal = 0, 
                           const int &in_coronal = 0, 
                           const int &in_axial = 0,
                           const int &n_in = 3 );

#endif