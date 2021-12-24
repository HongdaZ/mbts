#ifndef CALCDESCR_H
#define CALCDESCR_H

#include <vector>
#include <list>

using std::vector;
using std::list;

// Calculate the shape descriptors for the slices
void calcDescr( const int &size_region,
                int *ptr_seg_copy, const int &len, 
                vector<int> &region,
                const vector<list<vector<int>>> &slices, 
                vector<double> &solidity, 
                vector<double> &avg_spread, 
                vector<double> &avg_round,
                const int *ptr_aidx, const int *ptr_nidx );

#endif