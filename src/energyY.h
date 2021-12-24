#ifndef ENERGYY_H
#define ENERGYY_H

#include <vector>

using std::vector;

// calculate energy for tumor regions
double energyY( const vector<int> &region,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                const double &alphak,
                const double &betak,
                const double &a, const double &b );
// Single point energyY for healty voxel
double energyY( const int &curr_label,
                const int &curr_idx,
                const double &mu,
                const double &sigma2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta );
// Single point energyY for tumor voxel
// having tumor neighbor
double energyY( const int curr_label,
                const int curr_idx,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                const double &alphak,
                const double &betak,
                const double &a, const double &b );
// Single point energyY for tumor without tumor neighbor
// or outlier
double energyY( const int &curr_label,
                const int &curr_idx,
                const double &mu,
                const double &mk1,
                const double &sigma2,
                const double &lambda2,
                const int *ptr_seg,
                const double *ptr_intst,
                const double &alphak,
                const double &betak,
                const double &a, const double &b );
                            

#endif