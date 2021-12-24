#ifndef INITPARM_H
#define INITPARM_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// Initialize parameters
void initParm( const int &n_health,
               vector<int> &region, vector<double> &theta, 
               const bool &first_run, vector<double> &health_parm,
               vector<double> &tumor_parm,
               int *ptr_seg, const double *ptr_m,
               const double *ptr_nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               const list<list<int>> &tumor_regions,
               const double *ptr_a, const double *ptr_b, const int &len, 
               const int &maxit );
// Initialize parameters for t1ce, flair and t2 images
void initParmHealth( const int &n_health,
                     vector<int> &region, vector<double> &theta,
                     vector<double> &health_parm,
                     int *ptr_seg, const double *ptr_m,
                     const double *ptr_nu2, const double *ptr_intst, 
                     const double *ptr_lambda2, const int *ptr_nidx,
                     const double *ptr_nintst, const double *ptr_alpha,
                     const double *ptr_beta, const int &len, 
                     const int &maxit );
#endif