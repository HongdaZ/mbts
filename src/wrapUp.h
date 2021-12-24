#ifndef WRAPUP_H
#define WRAPUP_H

// Wrap up results for HGG
void wrapUp( const int &len, 
             int *ptr_hemorrhage,
             int *ptr_necrosis,
             int *ptr_enh,
             int *ptr_edema,
              int *ptr_seg,
              int *ptr_tumor );
// Wrap up results for last step
void wrapUp( const int &len, int *ptr_edema_regions,
             int *ptr_hemorrhage,
             int *ptr_necrosis,
             int *ptr_enh,
             int *ptr_edema,
             const int *ptr_flair,
             int *ptr_seg, int *ptr_tumor );

#endif