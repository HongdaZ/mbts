#ifndef GETPARM_H
#define GETPARM_H

#include <vector>

using std::vector;

// get existing parameters
void getParm( double &mu, double &sigma2, vector<double> &theta,
             const vector<double> &x_parm, const int &cidx );
void getParm( double &mu, double &sigma2,
              const vector<double> &outl_parm, const int &cidx );
#endif