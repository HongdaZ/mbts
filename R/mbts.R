#' A Three-stage Procedure for Multimodal Brain Tumor Segmentation
#' 
#' \code{mbts} analyzes T1ce, FLAIR and T2 images and segments the images 
#' into brain tissues.
#' @author Hongda Zhang
#' @useDynLib mbts
#' @importFrom stats kmeans predict smooth.spline quantile
#' @importFrom oro.nifti readNIfTI writeNIfTI nifti
#' @param patient A vector of file names of FLAIR, T1, T1ce and T2 images.
#' @param out The folder to save the intermediate and final results.
#' @param  infolder The folder which includes the multimodal images.
#' @param shrink Whether to shrink the original image by a factor of 2
#' in each dimension of 3D space.
#' @param delta A list of values related to unary potentials for 
#' segmenting T1ce, FLAIR and T2 images and for spliting non-enhancing 
#' tumor and edema.
#' @param delta_factor A factor for adjusting unary potentials. The larger the
#' values, the smaller the resulting unary potentials.
#' @param gamma A list of parameters for pairwise potentials. The larger the 
#' values the stronger the spatial correlation between neighboring voxels.
#' @param alpha,beta  Shape and scale parameter of inverse-gamma 
#' prior of variances.
#' @param lambda2 Prior variance of interaction parameters.
#' @param a Shape parameter of prior distribution of mean intensity of 
#' bright signals.
#' @param nu2 Variance of mean intensity of distinguishable classes.
#' @param maxit Maximum iterations to run for getting the segmentation results.
#' @param min_enh Minimum number of enhancing tumor voxels 
#' a high-glioma patient has.
#' @param min_enh_enc Minimum number of enhancing tumor voxels enclosed in
#' tumor regions a high-glioma patient has.
#' @param max_prop_enh_enc Maximal proportion of enhancing tumor voxels 
#' enclosed in tumor regions.
#' @param max_prop_enh_slice Maximal proportion of enhancing tumor voxels 
#' enclosed in a 2D slice of tumor regions.
#' @param min_tumor Minimal size of a tumor region.
#' @param spread_add The maximal form factor a region to be added can have.
#' @param spread_rm The minimal form factor a region to be removed can have.
#' @param trim1_spread,trim1_round Form factor and roundness for trimming the 
#' regions for the first time.
#' @param remove2d_spread,remove2d_round 2D form factor and roundness for 
#' removing 2D slices.
#' @param spread_trim,round_trim Form factor and roundness for trimming the 3D 
#' regions.
#' @param on_flair_prop,on_flair_hull_prop,on_flair_nt_prop Proportions for 
#' extending whole tumor regions in FLAIR images.
#' @param last_rm_solidity,last_rm_spread,last_rm_round Solidity, form factor
#' and roundness for removing the 3D tumor regions for the last time.
#' @param last_trim_spread,last_trim_round Form factor and roundness for
#' trimming the tumor regions for the last time.
#' @param last_trim_rm_spread,last_trim_rm_round Form factor and roundness
#' for removing tumor regions after trimming the tumor regions for the 
#' last time.
#' @param csf_check Whether to validate some of the necrosis 
#' voxels are actually CSF. 
#' @export 
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
                    list( t1ce = 40L, flair = 40L,
                          t2 = 40L, fthr = 40L ),
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
    fthr <- NULL
    fthr$delta <- delta$fthr
    fthr$gamma <- gamma$fthr
    fthr$alpha <- alpha$fthr
    fthr$beta <- beta$fthr
    fthr$lambda2 <- lambda2$fthr
    fthr$nu2 <- nu2$fthr
    fthr$maxit <- maxit$fthr
    
    post( patient, out, infolder, fthr, a,
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
