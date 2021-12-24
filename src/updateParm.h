#ifndef UPDATEPARM_H
#define UPDATEPARM_H

#include <vector>

#include "updateMu.h"
#include "updateTS.h"
#include "updateSigma.h"

using std::vector;

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const vector<int> &region,
                 const double &m,
                 const double &nu2,
                 const double *ptr_intst,
                 const int &curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double &alphal,
                 const double &betal,
                 const int &maxit );
// update parameters for tumor regions
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const vector<int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 const int &curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 const int &maxit );
// update parameters for outliers or single point tumor regions
void updateParm( double &mu, double &sigma2, 
                 const int idx,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 int *ptr_seg,
                 const double alphal,
                 const double betal,
                 const int &maxit );

#endif