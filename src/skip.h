#ifndef SKIP_H
#define SKIP_H

#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "nbrLabel.h"

using std::vector;

// whether skip current voxel
bool skip( const int &idx, const int *ptr_seg, const int *ptr_nidx );

#endif