#ifndef PRINTPARM_H
#define PRINTPARM_H

#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <map>

using std::vector;
using std::map;

void printParm( map<int, vector<double>> &the_parm );

#endif