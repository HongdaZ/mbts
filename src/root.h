#include <vector>

using std::vector;

// solve cubic equation x^3 + a*x^2 + b*x + c = 0
// x - array of size 3
// return 3: 3 real roots x[0], x[1], x[2]
// return 1: 1 real root x[0] and pair of complex roots: x[1]±i*x[2]

int root( vector<double> &x, double a, double b, double c );		
