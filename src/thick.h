#ifndef THICK_H
#define THICK_H

// get thickness of ptr_seg2 against ptr_seg1
int thick( const int *ptr_seg1, const int &label1, 
           int *ptr_seg2, const int &label2,
           const int &len, const int *ptr_nidx );
#endif