#include "expRatio.h"

double expRatio( const double &r0, const double &h ) {
  double ratio = 0;
  ratio = r0 * r0 / ( r0 * r0 + h * h );
  return ratio;
}