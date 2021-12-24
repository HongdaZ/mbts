#include "tissueType.h"
#include "wrapUp.h"

// Wrap up results for HGG
void wrapUp( const int &len, 
             int *ptr_hemorrhage,
             int *ptr_necrosis,
             int *ptr_enh,
             int *ptr_edema,
             int *ptr_seg,
             int *ptr_tumor ) {
  for( int i = 0; i < len; ++ i ) {
    if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
      ptr_necrosis[ 2 * i ] = 0;
      ptr_enh[ 2 * i ] = 0;
      ptr_edema[ 2 * i ] = 0;
    } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
      ptr_hemorrhage[ 2 * i ] = 0;
      ptr_enh[ 2 * i ] = 0;
      ptr_edema[ 2 * i ] = 0;
    } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ptr_seg[ 2 * i ] = Seg::SET;
      ptr_tumor[ 2 * i ] = 1;
      ptr_hemorrhage[ 2 * i ] = 0;
      ptr_necrosis[ 2 * i ] = 0;
      ptr_edema[ 2 * i ] = 0;
    } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
      ptr_seg[ 2 * i ] = Seg::SED;
      ptr_tumor[ 2 * i ] = 1;
      ptr_hemorrhage[ 2 * i ] = 0;
      ptr_necrosis[ 2 * i ] = 0;
      ptr_enh[ 2 * i ] = 0;
    }
  }
}

// Wrap up results for last step
void wrapUp( const int &len, int *ptr_edema_regions,
             int *ptr_hemorrhage,
             int *ptr_necrosis,
             int *ptr_enh,
             int *ptr_edema,
             const int *ptr_flair,
             int *ptr_seg, int *ptr_tumor ) {
  int remainder = 0;
  for( int i = 0; i < len; ++ i ) {
    remainder = ptr_edema_regions[ 2 * i ] % 10;
    if( remainder == 1 ) {
      if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_necrosis[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
        ptr_seg[ 2 * i ] = Seg::SET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_necrosis[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        
        if( ptr_flair[ 2 * i ] == Flair::FCSF ||
            ptr_flair[ 2 * i ] == Flair::FWM || 
            ptr_flair[ 2 * i ] == Flair::FGM ) {
          ptr_seg[ 2 * i ] = Seg::SNET;
          ptr_edema[ 2 * i ] = 0;
          ptr_necrosis[ 2 * i ] = Tumor::NCR;
        } else {
          ptr_seg[ 2 * i ] = Seg::SED;
        }
      }
    } else {
      if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_necrosis[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
        ptr_seg[ 2 * i ] = Seg::SET;
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_necrosis[ 2 * i ] = 0;
        ptr_edema[ 2 * i ] = 0;
      } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
        ptr_seg[ 2 * i ] = Seg::SED;
        ptr_tumor[ 2 * i ] = 1;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_necrosis[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
      }
    }
  }
}