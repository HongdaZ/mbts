#ifndef ASSIGNIDX_H
#define ASSIGNIDX_H

// Assign index of neighboring voxels
void assignIdx( const int &i, const int &j, const int &k,
                const int &nr, const int &nc, const int &ns,
                const int &n_nbr, const int &n_valid,
                const int *fidx, int *ptr_nidx );
#endif