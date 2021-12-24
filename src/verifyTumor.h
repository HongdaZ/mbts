#ifndef VERIFYTUMOR_H
#define VERIFYTUMOR_H

#include <map>
#include <list>

using std::map;
using std::list;

// check if tumor regions match seg
bool verifyTumor( map<int, list<int>> &tumor_regions, const int *ptr_seg,
                  const int len );

#endif