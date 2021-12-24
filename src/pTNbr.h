#ifndef PTNBR_H
#define PTNBR_H

#include <vector>

using std::vector;

double pTNbr( const vector<int> &region, const int *ptr_seg1, 
            const int &label1, const int *ptr_nidx );
double pTNbr2D( const vector<int> &region, const int *ptr_seg1, 
                const int &label1, const int *ptr_nidx, 
                const int *ptr_aidx, const int &plane );

#endif