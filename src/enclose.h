#ifndef ENCLOSE_H
#define ENCLOSE_H

#include <vector>
#include <list>

using std::vector;
using std::list;

// Find part of slices2 enclosed in slices1
void enclose( int *ptr_seg, const int &len,
              const vector<list<vector<int>>> &slices1,
              const vector<list<vector<int>>> &slices2,
              vector<double> &volume,
              const int &in_sagittal = 0, 
              const int &in_coronal = 0, 
              const int &in_axial = 0,
              const int &n_in = 3 );

#endif