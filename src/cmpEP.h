#ifndef CMPEP_H
#define CMPEP_H

#include <vector>
#include <map>
#include <list>

#include "updateParm.h"
#include "energyY.h"
#include "energyX.h"

using std::vector;
using std::map;
using std::list;

// compare energy for prediction
void cmpEP( const int &n_health,
            const list<int> &update_parm,
            vector<int> &region, const int &idx, const int &sc,
            const vector<int> &regions_whole,
            const vector<int> &regions_sub,
            vector<int> &tumor_labels, vector<int> &outl_labels,
            const vector<double> &health_parm,
            vector<double> &tumor_parm, 
            vector<double> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2,
            // outlier_parm( 3, 0 )
            vector<double> &outlier_parm,
            // theta( 6, 0 )
            vector<double> &theta,
            // tmp_parm( 9, 0 )
            vector<double> &tmp_parm,
            // out_theta( 6, 0 )
            vector<double> &out_theta,
            // new_out_parm( 2, 0 )
            vector<double> &new_out_parm,
            // whole_parm( 8, 0 )
            vector<double> &whole_parm,
            // label_whole_parm( 9, 0 )
            vector<double> &label_whole_parm,
            // used to store parameters of whole region and sub-regions
            vector<double> &region_parm,
            int &n_tumor, int &n_outl,
            const int &n_region,
            list<list<int>> &tumor_regions,
            const int &n_row,
            const list<int> &tumor_label );

#endif