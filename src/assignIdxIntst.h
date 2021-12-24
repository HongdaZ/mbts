#ifndef ASSIGNIDXINTST_H
#define ASSIGNIDXINTST_H

// Assign value to neighbor idx and neighbor intst
void assignIdxIntst( const int &i, const int &j, const int &k,
                     const int &nr, const int &nc, const int &ns,
                     const int &n_nbr, const int &n_valid,
                     const int *ptr_fidx, int *ptr_idx, double *ptr_intst,
                     const double *ptr_vintst );

#endif