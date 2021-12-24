#include <R.h>
#include <Rinternals.h>

// Assign value to neighbor idx and neighbor intst
void assignIdxIntst( const int &i, const int &j, const int &k,
                     const int &nr, const int &nc, const int &ns,
                     const int &n_nbr, const int &n_valid,
                     const int *ptr_fidx, int *ptr_idx, double *ptr_intst,
                     const double *ptr_vintst ) {
  int offset[ 18 ] = { 0, 0, -1, 
                       0, -1, 0,
                       -1, 0, 0,
                       1, 0, 0, 
                       0, 1, 0, 
                       0, 0, 1 };
  for( int l = 0; l < n_nbr; ++ l ) {
    ptr_idx[ n_nbr * n_valid + l ] = NA_INTEGER;
    ptr_intst[ n_nbr * n_valid + l ] = R_NaN;
  }
  int i_, j_, k_, tmp_index, tmp_vidx;
  for( int l = 0; l < n_nbr; ++ l ) {
    i_ = i + offset[ 3 * l + 0 ];
    j_ = j + offset[ 3 * l + 1 ];
    k_ = k + offset[ 3 * l + 2 ];
    if( i_ < nr && i_ >= 0 && 
        j_ < nc && j_ >= 0 &&
        k_ < ns && k_ >= 0 ) {
      tmp_index = nc * nr * k_ + nr * j_ + i_ + 1;
      
      if( ptr_fidx[ tmp_index - 1 ] != NA_INTEGER ) {
        tmp_vidx = ptr_fidx[ tmp_index - 1 ];
        ptr_idx[ n_nbr * n_valid + l ] = tmp_vidx;
        ptr_intst[ n_nbr * n_valid + l ] =
          ptr_vintst[ tmp_vidx - 1 ];
      }
    }
  }
  return;
}