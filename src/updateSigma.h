#ifndef UPDATESIGMA_H
#define UPDATESIGMA_H

// update sigma for outliers
void updateSigma( const int idx,
                  const double mu,
                  const double *ptr_intst,
                  const double alphal,
                  const double betal,
                  double &sigma2 );

#endif