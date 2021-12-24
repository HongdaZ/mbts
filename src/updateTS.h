#ifndef UPDATETS_H
#define UPDATETS_H

#include <list>
#include <vector>

using std::list;
using std::vector;

// region starts from 1
// updateTheta and sigma2 for health and tumorous regions
void updateTS( const int &nrow,
               const double &mu,
               double &sigma2,
               const double &lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               vector<double> &theta,
               const double &alphal,
               const double &betal,
               double *yl, double *yln, double *ylna,
               const double *yln_,
               const double *yln_i,
               const double *yl_ );
                            

#endif