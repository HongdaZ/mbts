#include "printParm.h"

void printParm( map<int, vector<double>> &the_parm ) {
  for( map<int, vector<double>>::iterator it = the_parm.begin();
       it != the_parm.end(); ++ it ) {
    int label = it->first;
    Rprintf( "label = %d; parm =", label );
    vector<double> &parm = it->second;
    for( vector<double>::iterator it_parm = parm.begin();
         it_parm != parm.end(); ++ it_parm ) {
      Rprintf( "%f, ", *it_parm );
    }
    Rprintf( "\n" );
  }
  return;
}
