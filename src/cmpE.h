#ifndef CMPE_H
#define CMPE_H

#include <vector>
#include <map>

#include "energyY.h"
#include "energyX.h"

using std::vector;
using std::map;

// compare energy for t1ce, flair or t2 images in estimation
void cmpE( const int &n_health,
           int idx, vector<double> &health_parm,
           int *ptr_seg, const int *ptr_nidx,
           const double *ptr_intst, const double *ptr_nintst,
           const double *ptr_delta, const double *ptr_gamma,
           vector<double> &theta );
#endif