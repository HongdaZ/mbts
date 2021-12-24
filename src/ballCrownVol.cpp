#include "ballCrownVol.h"

double ballCrownVol( const double &r0, const double &h ) {
  const double PI = 3.14159265358979323846;
  double v;
  v = PI * h * ( 3 * r0 * r0 + h * h ) / 6;
  return v;
}