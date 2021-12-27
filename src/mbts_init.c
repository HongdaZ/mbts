#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "header.h"

// Declare functions
SEXP indexMat( SEXP img, SEXP label );
SEXP infoMat( SEXP img, SEXP label );
SEXP changeB( SEXP label, SEXP image1, SEXP q1, SEXP image2, 
              SEXP q2, SEXP k  );
SEXP postProcess( SEXP post_data, SEXP min_enh,
                  SEXP min_enh_enc,
                  SEXP max_prop_enh_enc,
                  SEXP max_prop_enh_slice,
                  SEXP min_tumor, SEXP spread_add,
                  SEXP spread_rm, 
                  SEXP trim1_spread, SEXP trim1_round,
                  SEXP remove2d_spread, SEXP remove2d_round,
                  SEXP spread_trim, SEXP round_trim,
                  SEXP on_flair_prop, 
                  SEXP on_flair_hull_prop,
                  SEXP on_flair_nt_prop,
                  SEXP last_rm_solidity, 
                  SEXP last_rm_spread, 
                  SEXP last_rm_round,
                  SEXP last_trim_spread,
                  SEXP last_trim_round,
                  SEXP last_trim_rm_spread,
                  SEXP last_trim_rm_round,
                  SEXP csf_check );
SEXP initPost( SEXP t1ce_image, SEXP flair_image, 
               SEXP t2_image, SEXP t1ce_intst );
SEXP changeA( SEXP img, SEXP label );
SEXP convexHull( SEXP points );
SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
              SEXP alpha, SEXP beta, SEXP lambda2, 
              SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit );
SEXP changeD( SEXP old, SEXP a, SEXP b );
SEXP inSide( SEXP p, SEXP poly );
SEXP priorMode( SEXP model );
SEXP est( SEXP model, SEXP delta, SEXP gamma, 
          SEXP alpha, SEXP beta, SEXP lambda2, 
          SEXP m, SEXP nu2, SEXP maxit );
SEXP changeC( SEXP old, SEXP image, SEXP q );
SEXP call_interface( SEXP x );
SEXP pred( SEXP model, SEXP delta, SEXP gamma,
           SEXP alpha, SEXP beta, SEXP lambda2,
           SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit );
SEXP estF( SEXP model, SEXP delta, SEXP gamma, 
           SEXP alpha, SEXP beta, SEXP lambda2, 
           SEXP m, SEXP nu2, SEXP maxit );
// Include all .c functions
static R_CallMethodDef callMethods[] = {
  { "indexMat", (DL_FUNC) &indexMat, 2 },
  { "infoMat", (DL_FUNC) &infoMat, 2 },
  { "changeB", (DL_FUNC) &changeB, 6 },
  { "postProcess", (DL_FUNC) &postProcess, 25 },
  { "initPost", (DL_FUNC) &initPost, 4 },
  { "changeA", (DL_FUNC) &changeA, 2 },
  { "convexHull", (DL_FUNC) &convexHull, 1 },
  { "estParm", (DL_FUNC) &estParm, 11 },
  { "changeD", (DL_FUNC) &changeD, 3 },
  { "inSide", (DL_FUNC) &inSide, 2 },
  { "priorMode", (DL_FUNC) &priorMode, 1 },
  { "est", (DL_FUNC) &est, 9 },
  { "changeC", (DL_FUNC) &changeC, 3 },
  { "call_interface", (DL_FUNC) &call_interface, 1 },
  { "pred", (DL_FUNC) &pred, 11 },
  { "estF", (DL_FUNC) &estF, 9 },
  { NULL, NULL, 0 }
};
/* This is called by the dynamic loader to register the routine. */
void R_init_mbts(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}