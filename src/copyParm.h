#ifndef COPYPARM_H
#define COPYPARM_H

#include <list>
#include <vector>

using std::list;
using std::vector;

// copy estimated parameters to SEXP
void copyParm( const int &n_health,
               vector<double> &health_parm,
               vector<double> &tumor_parm,
               vector<double> &outl_parm,
               double *ptr_res_parm, const int &nrow, 
               const list<list<int>> &tumor_regions, 
               const vector<int> &outl_labels, const int &len );
void copyParmHealth( const int &n_health, vector<double> &health_parm,
                     double *ptr_res_parm, const int &nrow );
#endif