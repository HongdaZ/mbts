#ifndef PERIMETER_H
#define PERIMETER_H

// Find perimeter of a 2D region
int perimeter( const int *ptr_seg, const int &label, 
               const int &len, const int *ptr_nidx,
               const int &plane );
// 3D version of the function above
int perimeter( const int *ptr_seg, const int &label,
               const int &len, const int *ptr_nidx );

#endif