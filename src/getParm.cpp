#include "getParm.h"

void getParm( double &mu, double &sigma2, vector<double> &theta,
             const vector<double> &x_parm, const int &cidx ) {
  mu = x_parm[ 8 * cidx ];
  sigma2 = x_parm[ 8 * cidx + 1 ];
  for( int j = 0; j < 6; ++ j ) {
    theta[ j ] = x_parm[ 8 * cidx + 2 + j ];
  }
  return;
}

void getParm( double &mu, double &sigma2,
              const vector<double> &outl_parm, const int &cidx ) {
  mu = outl_parm[ 2 * cidx ];
  sigma2 = outl_parm[ 2 * cidx + 1 ];
  return;
}