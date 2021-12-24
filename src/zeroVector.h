#ifndef ZEROVECTOR_H
#define ZEROVECTOR_H

#include <vector>

using std::vector;

// assign zero to all elements of the vector
template< typename T >
void zeroVector( vector<T> &v ) {
  for( typename vector<T>::iterator it = v.begin();
       it != v.end(); ++ it ) {
    *it = 0;
  }
  return;
}
template< typename T >
void zeroVector( T *ptr, const int &len ) {
  const int length = 2 * len;
  for( int i = 0; i < length; ++ i ) {
    ptr[ i ] = 0;
  }
}
#endif