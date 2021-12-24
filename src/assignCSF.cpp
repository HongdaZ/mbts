#include "assignCSF.h"
#include "tissueType.h"

// Find CSF
void assignCSF( int *ptr_csf, const int *ptr_t1ce,
                const int *ptr_flair, const int *ptr_tumor, 
                const int &len ) {
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_t1ce[ 2 * i ] == T1ce::T1CSF ||
        ptr_flair[ 2 * i ] == Flair::FCSF ) &&
        ptr_tumor[ 2 * i ] == 0 ) {
      ptr_csf[ 2 * i ] = 1;
    } else {
      ptr_csf[ 2 * i ] = 0;
    }
  }
}