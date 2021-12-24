#ifndef UPDATEBETA_H
#define UPDATEBETA_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// update beta for healthy and tumor regions
void updateBeta( const int &n_health,
                 double *ptr_beta, const double *ptr_alpha, 
                 const vector<double> &health_parm,  
                 const list<list<int>> &tumor_regions,
                 const vector<double> &tumor_parm );
// update beta for t1ce, flair or t2 images
void updateBeta( const int &n_health,
                 double *ptr_beta, const double *ptr_alpha, 
                 const vector<double> &health_parm );
#endif