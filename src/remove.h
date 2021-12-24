#ifndef REMOVE_H
#define REMOVE_H

#include <vector>
#include <limits>

using std::vector;

// Remove extra regions of tumor
void remove( vector<int> &region, 
             int *ptr_seg2_copy,
             int *ptr_seg_copy,
             const int *ptr_aidx, 
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis, 
             int *ptr_enh, int *ptr_edema, const int &m_tumor, 
             const int &m_enh, const int &m_enh_enc,
             const double &spread_factor,
             const int len, const int &nr, const int &nc, const int &ns );
// Remove small regions (calls removeSmall)
void remove( vector<int> &region, 
             int *ptr_seg2_copy,
             int *ptr_seg_copy,
             const int *ptr_aidx, 
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis, 
             int *ptr_enh, int *ptr_edema, const int &m_tumor, 
             const int &m_enh, const int &m_enh_enc, const int len, 
             const int &nr, const int &nc, const int &ns,
             const double &rm_spread = 2, 
             const double &rm_round = 
               std::numeric_limits<double>::infinity() );
#endif