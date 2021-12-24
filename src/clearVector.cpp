#include "clearVector.h"

// erase all elements from a container
template <typename T >
void clearVector( T &v );
template <>
void clearVector<list<int>>( list<int> &v ){
  v.clear();
  return;
};
template <>
void clearVector<list<double>>( list<double> &v ){
  v.clear();
  return;
};
template <>
void clearVector<vector<int>>( vector<int> &v ){
  while( v.size() > 0 ) {
    v.pop_back();
  }
  return;
};
template <>
void clearVector<vector<double>>( vector<double> &v ){
  while( v.size() > 0 ) {
    v.pop_back();
  }
  return;
};