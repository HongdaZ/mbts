#ifndef EXCLUDENECROSIS_H
#define EXCLUDENECROSIS_H

#include <vector>

using std::vector;
// Exclude necrosis regions with mean t1ce intensity < mean t1ce csf
void excludeNecrosis( const int &len, const int *ptr_nidx,
                      int *ptr_new, 
                      const int *ptr_t1ce_csf,
                      const double *ptr_t1ce_intst,
                      int *ptr_add_necrosis,
                      vector<int> &add_necrosis_region,
                      const double &mean_csf );

#endif