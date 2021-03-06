#' Postprocessing Stage of the Procedure
#' 
#' \code{post} postprocesses the images.
#' @author Hongda Zhang
#' @param patient A vector of file names of FLAIR, T1, T1ce and T2 images.
#' @param out The folder to save the intermediate and final results.
#' @param infolder The folder which includes the multimodal images.
#' in each dimension of 3D space.
#' @param fthr Parameters for further splitting ED and NET.
#' @param a Shape parameter of prior distribution of mean intensity of 
#' bright signals.
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
# postprocessing after segmentation
post <- function( patient, out = "SEG", infolder = "N4ITK433Z",
          fthr = list( delta = c( 0, 0, 4, 5 ),
                       gamma = 0.8,
                       alpha = c( 10, 10 ),
                       beta = c( 1, 1 ),
                       lambda2 = c( 1, 1 ),
                       nu2 = rep( .25, 2 ),
                       maxit = 40L ),
          a = 2,
          min_enh = 500L,
          min_enh_enc = 2000L,
          max_prop_enh_enc = .1,
          max_prop_enh_slice = .2,
          min_tumor = 20000L,
          spread_add = 8,
          spread_rm = 9,
          trim1_spread = 8,
          trim1_round = 17,
          remove2d_spread = 25,
          remove2d_round = 25,
          spread_trim = 7,
          round_trim = 16,
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
  if( is.null( last_trim_spread ) ) {
    last_trim_spread <- spread_trim
  }
  if( is.null( last_trim_round ) ) {
    last_trim_round <- round_trim
  }
  ## Read segmentation results
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  out_t1ce_seg <- gsub( "_flair.nii.gz", "_t1ce_seg", outfile )
  t1ce_image <- readNIfTI( out_t1ce_seg, reorient = FALSE )@.Data
  out_t1ce_intst <- gsub( "_flair.nii.gz", "_t1ce_norm", outfile )
  t1ce_intst <- readNIfTI( out_t1ce_intst, reorient = FALSE )@.Data
  out_flair_seg <- gsub( "_flair.nii.gz", "_flair_seg", outfile )
  flair_image <- readNIfTI( out_flair_seg, reorient = FALSE )@.Data
  out_t2_seg <- gsub( "_flair.nii.gz", "_t2_seg", outfile )
  t2_image <- readNIfTI( out_t2_seg, reorient = FALSE )@.Data
  ## Initialize data for postprocessing
  post_data <- initPost( t1ce_image, flair_image, t2_image,
                         t1ce_intst )
  # sink( '/media/hzhang/ZHD-P1/result/output.txt' )
  post_seg <- postProcess( post_data, min_enh, min_enh_enc,
                           max_prop_enh_enc, max_prop_enh_slice,
                           min_tumor, spread_add, spread_rm,
                           trim1_spread, trim1_round, remove2d_spread,
                           remove2d_round, spread_trim, round_trim,
                           on_flair_prop, on_flair_hull_prop,
                           on_flair_nt_prop,
                           last_rm_solidity, 
                           last_rm_spread,
                           last_rm_round,
                           last_trim_spread, last_trim_round, 
                           last_trim_rm_spread,
                           last_trim_rm_round,
                           csf_check )
  # sink()
  if( sum( post_seg$image == 2 |
           post_seg$image == 1 |
           post_seg$image == 4 , na.rm = T ) < 2000 ) {
    ## Needs to segment t2 again
    post_seg$image[ post_seg$image == 6 ] <- NA_integer_
    post_seg$image[ post_seg$image == 5 ] <- NA_integer_
    out_t2_data <- gsub( "_flair.nii.gz", "_t2.RData", outfile )
    load( out_t2_data )
    t2_shift <- NULL
    new_delta_t2 <- t2_shift
    new_delta_t2[ 3 ] <- new_delta_t2[ 3 ] - 0.5
  } else {
    new_delta_t2 <- NULL
    ## Seg csf
    if( length( post_seg$csf_code ) > 0 ) {
      for( i in 1 : length( post_seg$csf_code ) ) {
        csf_code <- post_seg$csf_code[ i ]
        if( sum( post_seg$csf == csf_code, na.rm = T ) > 10 ) {
          ## With CSF inside tumor
          out_t1ce_norm <- gsub( "_flair.nii.gz", "_t1ce_norm", outfile )
          t1ce_intst <- readNIfTI( out_t1ce_norm, 
                                   reorient = FALSE )@.Data
          ## Furtherly segment CSF 
          sub_csf_image <- array( NA_integer_, dim = dim( t1ce_intst ) )
          sub_csf_idx <- post_seg$image == 1 | 
            post_seg$image == 5 |
            post_seg$csf == csf_code
          sub_csf_idx[ is.na( sub_csf_idx ) ] <- FALSE
          sub_csf_image[ sub_csf_idx ] <- post_seg$image[ sub_csf_idx ]
          further_data <- splitFthrC( sub_csf_image, t1ce_intst )
          m <- further_data$m
          further_model <-initFther( further_data$label, 
                                     further_data$intst )
          further_seg <- est( further_model, fthr$delta[ 1 : 2 ], 
                              fthr$gamma[ 1 : 2 ],
                              fthr$alpha[ 1 : 2 ], 
                              fthr$beta[ 1 : 2  ], 
                              fthr$lambda2[ 1 : 2 ],
                              m, fthr$nu2[ 1 : 2 ], 
                              fthr$maxit )
          m <- further_seg$parm[ 2, ]
          sigma2 <- further_seg$parm[ 3, ]
          if( ( m[ 2 ] - m[ 1 ] ) /  fthr$delta[ 3 ] > 
              sqrt( min( sigma2 ) ) ) {
            necrosis_idx <- sub_csf_image == 6 &
              further_seg$image == -2
            necrosis_idx[ is.na( necrosis_idx ) ] <- FALSE
            post_seg$image[ necrosis_idx ] <- 1L
          }
        }
      }
    }
     
    post_seg$image[ post_seg$image == 6 ] <- NA_integer_
    post_seg$image[ post_seg$image == 5 ] <- NA_integer_
    # sink()
    if( length( post_seg$edema_code ) > 0 ) {
      for( i in 1 : length( post_seg$edema_code ) ) {
        edema_code <- post_seg$edema_code[ i ]
        remainder <- edema_code %% 10
        if( remainder != 3 ) {
          if( sum( post_seg$edema == edema_code, na.rm = T) > 10 ) {
            out_t2_norm <- gsub( "_flair.nii.gz", "_t2_norm", outfile )
            t2_intst <- readNIfTI( out_t2_norm, reorient = FALSE )@.Data
            ## Furtherly segment edema
            sub_edema_image <- array( NA_integer_, 
                                    dim = dim( t2_intst ) )
            sub_edema_idx <- post_seg$edema == edema_code
            sub_edema_idx[ is.na( sub_edema_idx ) ] <- FALSE
            sub_edema_image[ sub_edema_idx ] <- 
              post_seg$image[ sub_edema_idx ]
            if( sum( sub_edema_image == 2 |
                     sub_edema_image == 4, na.rm = T) > 10) {
              further_data <- splitFthrE( sub_edema_image, t2_intst )
              m <- further_data$m
              further_model <-initFther( further_data$label, 
                                         further_data$intst )
              further_seg <- estF( further_model, fthr$delta, fthr$gamma,
                                   fthr$alpha, fthr$beta, fthr$lambda2,
                                   m, fthr$nu2, fthr$maxit )
              m <- further_seg$parm[ 2, ]
              sigma2 <- further_seg$parm[ 3, ]
              if( ( m[ 2 ] - m[ 1 ] ) /  fthr$delta[ 4 ] > 
                  sqrt( min( sigma2 ) ) ) {
                edema_idx <- sub_edema_image == 2
                edema_idx[ is.na( edema_idx ) ] <- FALSE
                post_seg$image[ edema_idx ] <- 
                  further_seg$image[ edema_idx ]
              }
            }
          }
        }
      }
    }
  }
  ## Export the results to .nii images
  post_seg$image[ is.na( post_seg$image ) ] <- 0
  out_post_seg <- gsub( "_flair.nii.gz", "_post_seg", outfile )
  writeNIfTI( nifti( post_seg$image, datatype = 2 ),
              filename = out_post_seg, gzipped = TRUE )
  out_new_delta_t2 <- gsub( "_flair.nii.gz", "_post.rds", outfile )
  saveRDS( new_delta_t2, out_new_delta_t2 )
}