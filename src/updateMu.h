#ifndef UPDATEMU_H
#define UPDATEMU_H

#include <list>
#include <vector>
#include <map>

using std::list;
using std::vector;
using std::map;

// update mu for healthy cells
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double nu2,
                 const double sum_theta,
                 const double sum_y );
// update mu for tumor cells
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double sum_theta,
                 const double sum_y );
// update mu for outliers
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double sum_y );

#endif