#ifndef SCPRED_H
#define SCPRED_H

#include <vector>
#include <list>

using std::vector;
using std::list; 

// split or combine for prediction
int scPred( int &n_region, list<int> &update_parm, 
            list<list<int>> &tumor_regions, 
            vector<int> &region, 
            const vector<int> &tumor_labels,
            int *ptr_label, const int *ptr_nidx, 
            const int &len, const int &start, 
            vector<int> &regions_whole, vector<int> &regions_sub,
            list<int> &tumor_nbr, list<int> &tumor_label );

#endif