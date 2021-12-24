#ifndef CLEARVECTOR_H
#define CLEARVECTOR_H

#include <vector>
#include <list>

using std::vector;
using std::list;

// erase all elements from a vector
template <typename T >
void clearVector( T &v );
template <>
void clearVector<list<int>>( list<int> &v );
template <>
void clearVector<list<double>>( list<double> &v );
template <>
void clearVector<vector<int>>( vector<int> &v );
template <>
void clearVector<vector<double>>( vector<double> &v );

#endif