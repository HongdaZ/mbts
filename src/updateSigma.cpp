#include <R.h>
#include <Rinternals.h>

#include "updateSigma.h"

void updateSigma( const int idx,
                  const double mu,
                  const double *ptr_intst,
                  const double alphal,
                  const double betal,
                  double &sigma2 ) {
  double intst = ptr_intst[ idx - 1 ];
  sigma2 = ( betal + pow( intst - mu, 2 ) / 2 ) / 
    ( 1 / (double)2 + alphal + 1 );
  return;
  
}
