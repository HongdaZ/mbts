#ifndef TISSUETYPE_H
#define TISSUETYPE_H

enum Tumor { HMG = 5, NCR = 6, NET = 7, ET = 4, ED = 2 };
enum T1ce { T1CSF = 1, T1GM = 2, T1WM = 3, T1TM = 4 };
enum Flair { FCSF = 1, FWM = 2, FGM = 3, FTM = 4 };
enum T2 { T2WM = 1, T2GM = 2, T2CSF = 4 };
enum Seg { SNET = 1, SET = 4, SED = 2, SCSF = 5, SECSF = 6 };
enum Plane { Sagittal = 0, Coronal = 1, Axial = 2 };

#endif