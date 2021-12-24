# postprocessing after segmentation
post <- function( patient, out = "SEG", infolder = "N4ITK433Z",
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
          further_seg <- est( further_model, delta$fthr[ 1 : 2 ], 
                              gamma$fthr[ 1 : 2 ],
                              alpha$fthr[ 1 : 2 ], 
                              beta$fthr[ 1 : 2  ], 
                              lambda2$fthr[ 1 : 2 ],
                              m, nu2$fthr[ 1 : 2 ], 
                              maxit$fthr )
          m <- further_seg$parm[ 2, ]
          sigma2 <- further_seg$parm[ 3, ]
          if( ( m[ 2 ] - m[ 1 ] ) /  delta$fthr[ 3 ] > 
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
              further_seg <- estF( further_model, delta$fthr, gamma$fthr,
                                   alpha$fthr, beta$fthr, lambda2$fthr,
                                   m, nu2$fthr, maxit$fthr )
              m <- further_seg$parm[ 2, ]
              sigma2 <- further_seg$parm[ 3, ]
              if( ( m[ 2 ] - m[ 1 ] ) /  delta$fthr[ 4 ] > 
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