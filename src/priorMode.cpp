#include <R.h>
#include <Rinternals.h>
#include "helper.h"

extern "C" SEXP priorMode( SEXP model );

SEXP priorMode( SEXP model ) {
  
  SEXP info = getListElement( model, "info" );
  SEXP intst = getListElement( info, "intst" );
  SEXP seg = getListElement( model, "seg" );
  
  const double *ptr_intst = REAL( intst );
  const int *ptr_seg = INTEGER( seg );
  int len = length( intst );
  
  double sum[ 4 ] = { 0 };
  int n[ 4 ] = { 0 };
  for( int i = 0; i < len; ++ i ) {
    switch( ptr_seg[ 2 * i ] ) {
    case -1 :
      n[ 0 ] += 1;
      sum[ 0 ] += ptr_intst[ i ];
      break;
    case -2 :
      n[ 1 ] += 1;
      sum[ 1 ] += ptr_intst[ i ];
      break;
    case -3 :
      n[ 2 ] +=1;
      sum[ 2 ] += ptr_intst[ i ];
      break;
    case 0 :
      break;
    default :
      n[ 3 ] += 1;
    sum[ 3 ] += ptr_intst[ i ];
    }
  }

    
  SEXP ans = PROTECT( allocVector( REALSXP, 4 ) );
  double *ptr_ans = REAL( ans );
  for( int i = 0; i < 4; ++ i ) {
    ptr_ans[ i ] = sum[ i ] / n[ i ];
  }
  UNPROTECT( 1 );
  return ans;
}
