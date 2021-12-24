#ifndef ASSIGNPARM_H
#define ASSIGNPARM_H

#include <vector>

using std::vector;

void assignParm( vector<double> &x_parm, const int &curr_label,
                  const vector<double> &parm );
void assignParm( vector<double> &x_parm, const int &curr_label, 
                  const double &mu, const double &sigma2,
                  const vector<double> &theta );

#endif