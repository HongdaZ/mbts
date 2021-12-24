# Markov random field model for brain tumor segmentation
mbts <- function( patient, out = "SEG", infolder = "N4ITK433Z", 
                  shrink = FALSE,
                  ## Always four numbers for delta
                  delta = 
                    list( t1ce = c( -1, 0, 8, 4 ),
                          flair = c( -0.45, 0, NA_real_, 4 ),
                          t2 = c( 2.65, 0, NA_real_, 4 ),
                          fthr = c( 0, 0, 4, 5 ) ),
                  delta_factor = 
                    list( t1ce = 1.75,
                          flair = 2.60,
                          t2 = 4.70 ),
                  gamma = list( t1ce = 0.8,
                                flair = 0.4,
                                t2 = 0.4,
                                fthr = 0.8 ),
                  ## #of healthy tissue types controlled by alpha
                  alpha = list( t1ce = c( 10, 10, 10, 10 ),
                                flair = c( 10, 10, 10, 10 ),
                                t2 = c( 10, 10, 10 ),
                                fthr = c( 10, 10 ) ),
                  beta = list( t1ce = c( 1, 1, 1, 1 ),
                               flair = c( 1, 1, 1, 1 ),
                               t2 = c( 1, 1, 1 ),
                               fthr = c( 1, 1 ) ),
                  lambda2 = list( t1ce = c( 1, 1, 1, 1 ),
                                  flair = c( 1, 1, 1, 1 ),
                                  t2 = c( 1, 1, 1 ),
                                  fthr = c( 1, 1 ) ),
                  a = 2,
                  nu2 = list( t1ce = rep( .25, 3 ),
                              flair = rep( .25, 3 ),
                              t2 = c( .25, .25 ),
                              fthr = rep( .25, 2 ) ),
                  maxit = 
                    list( t1ce = 10L, flair = 1L,
                          t2 = 1L, fthr = 40L ),
                  min_enh = 500L,
                  min_enh_enc = 2000L,
                  max_prop_enh_enc = .1,
                  max_prop_enh_slice = .2,
                  min_tumor = 2000L,
                  spread_add = 10,
                  spread_rm = 9,
                  trim1_spread = 10,
                  trim1_round = 18,
                  remove2d_spread = 25,
                  remove2d_round = 25,
                  spread_trim = 14,
                  round_trim = 20,
                  on_flair_prop = 1.5,
                  on_flair_hull_prop = 0.3,
                  on_flair_nt_prop = 0.3,
                  last_rm_solidity = 2,
                  last_rm_spread = 16,
                  last_rm_round = 16,
                  last_trim_spread = NULL,
                  last_trim_round = NULL,
                  last_trim_rm_spread = 2,
                  last_trim_rm_round = 10000,
                  csf_check = 0L ) {
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  out_new_delta_t2 <- gsub( "_flair.nii.gz", "_post.rds", outfile )
  redo <- TRUE
  if( file.exists( out_new_delta_t2 ) ) {
    new_delta_t2 <- readRDS( out_new_delta_t2 )
    if( is.null( new_delta_t2 ) ) {
      redo <- FALSE
    } else {
      if( new_delta_t2[ 3 ] < 2 ) {
        redo <- FALSE
      } else {
        delta$t2 <- new_delta_t2
      }
    }
  }
  while( redo ) {
    segment( patient, out, infolder, delta, delta_factor,
             gamma, alpha,
             beta, lambda2, a, nu2, maxit, redo, shrink )
    post( patient, out, infolder, delta, delta_factor,
          gamma, alpha,
          beta, lambda2, a, nu2, maxit, 
          min_enh, min_enh_enc, max_prop_enh_enc, 
          max_prop_enh_slice,
          min_tumor, spread_add, spread_rm,
          trim1_spread, trim1_round, remove2d_spread,
          remove2d_round, spread_trim, round_trim, 
          on_flair_prop, on_flair_hull_prop, on_flair_nt_prop,
          last_rm_solidity, last_rm_spread, last_rm_round,
          last_trim_spread, last_trim_round, last_trim_rm_spread,
          last_trim_rm_round, csf_check )
    if( file.exists( out_new_delta_t2 ) ) {
      new_delta_t2 <- readRDS( out_new_delta_t2 )
      if( is.null( new_delta_t2 ) ) {
        redo <- FALSE
      } else {
        if( new_delta_t2[ 3 ] < 2 ) {
          redo <- FALSE
        } else {
          delta$t2 <- new_delta_t2
        }
      }
    }
    break
  }
}
