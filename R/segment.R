# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, out = "SEG", infolder = "N4ITK433Z",
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
                     redo = TRUE,
                     shrink = FALSE
                     ) {
  images <- readImage( patient )
  if( shrink ) {
    images <- shrink( images )
  }
  ## Output files
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  
  out_t1ce_seg <- gsub( "_flair.nii.gz", "_t1ce_seg.nii.gz", outfile )
  out_t1ce_norm <- gsub( "_flair.nii.gz", "_t1ce_norm.nii.gz", outfile )
  out_t1ce_data <- gsub( "_flair.nii.gz", "_t1ce.RData", outfile )
  if( file.exists( out_t1ce_seg ) & 
      file.exists( out_t1ce_norm ) &
      file.exists( out_t1ce_data ) ) {
    t1ce_image <-  readNIfTI( out_t1ce_seg, reorient = FALSE )@.Data
  } else {
    ## split t1ce images to CSF & necrosis, grey matter and white matter
    t1ce_data <- splitT1ce3( images$t1ce, images$flair )
    m <- t1ce_data$m
    t1ce_model <- initEst( t1ce_data$label, t1ce_data$t1ce )
    ## estimate parameters of t2 images without tumor and CSF & necrosis
    t1ce_seg <- est( t1ce_model, delta$t1ce, gamma$t1ce,
                     alpha$t1ce[ 1 : 3 ], beta$t1ce[ 1 : 3 ], 
                     lambda2$t1ce[ 1 : 3 ], 
                     m, nu2$t1ce, 40L )
    t1ce_res <- t1ce_seg$parm[ c( 2 : 9 ), ]
    ## Update delta[ 1, 2 ] for t1ce
    shift_t1ce1 <- delta$t1ce[ 3 ]
    shift_t1ce2 <- updateDelta12( t1ce_res[ 1, 3 ], 
                                  t1ce_res[ 2, 3 ],
                                  t1ce_res[ 1, 1], 
                                  t1ce_res[ 2, 1 ], 
                                  shift_t1ce1 )
    delta$t1ce[ c( 1, 2 ) ] <- delta$t1ce[ c( 1, 2 ) ] - 
      shift_t1ce2
    t1ce_shift <- c( delta$t1ce[ 1 : 2 ], shift_t1ce1, delta$t1ce[ 4 ] )
    delta$t1ce[ 3 ] <- delta$t1ce[ 2 ] + 
      ( t1ce_shift[ 3 ] / delta_factor$t1ce ) ^ 2 / 2
    
    ## Export normalized images
    t1ce_intst <- t1ce_data$t1ce
    out_t1ce_norm <- gsub( ".nii.gz", "", out_t1ce_norm )
    writeNIfTI( nifti( t1ce_intst, datatype = 64 ) ,
                out_t1ce_norm,
                gzipped = T )
    save( t1ce_res, t1ce_shift, file = out_t1ce_data )
    ## update beta
    sigma2 <- t1ce_seg$parm[ 3, ]
    beta$t1ce[ 1 : 3 ] <- ( alpha$t1ce[ 1 : 3 ] + 1 ) * ( sigma2 )
    t1ce_data <- split4( t1ce_data$t1ce, t1ce_seg, t1ce_shift[ 3 ] )
    ## update m
    m <- t1ce_data$m
    b <- getB( m, a )
    t1ce_model <- initEst( t1ce_data$label, t1ce_data$intst )
    t1ce_seg <- pred( t1ce_model, delta$t1ce, 
                      gamma$t1ce, alpha$t1ce, beta$t1ce, 
                      lambda2$t1ce, a, b, m, nu2$t1ce, maxit$t1ce )
    ## t1ce
    t1ce_image <- t1ce_seg$image
    t1ce_image[ is.na( t1ce_image ) ] <- 0L
    t1ce_image[ t1ce_image >= 1 ] <- 4L
    t1ce_image[ t1ce_image <= -4 ] <- 4L
    t1ce_image[ t1ce_image == -1 ] <- 1L
    t1ce_image[ t1ce_image == -2 ] <- 2L
    t1ce_image[ t1ce_image == -3 ] <- 3L
    ## Get the initial results
    out_t1ce_seg <- gsub( ".nii.gz", "", out_t1ce_seg )
    writeNIfTI( nifti( t1ce_image, datatype = 2 ),
                filename = out_t1ce_seg, gzipped = TRUE )
  }
  
  out_flair_seg <- gsub( "_flair.nii.gz", "_flair_seg.nii.gz", outfile )
  out_flair_norm <- gsub( "_flair.nii.gz", "_flair_norm.nii.gz", outfile )
  out_flair_data <- gsub( "_flair.nii.gz", "_flair.RData", outfile )
  if( file.exists( out_flair_seg ) & 
      file.exists( out_flair_norm ) &
      file.exists( out_flair_data ) ) {
    flair_image <-  readNIfTI( out_flair_seg, reorient = FALSE )@.Data
  } else {
    ## split flair images to CSF & necrosis, white matter 
    ## and grey matter
    flair_data <- splitFlair3( images$flair, t1ce_image )
    m <- flair_data$m
    flair_model <- initEst( flair_data$label, flair_data$flair )
    ## estimate parameters of t1ce or FLAIR images without tumor
    flair_seg <- est( flair_model, delta$flair, gamma$flair, 
                      alpha$flair[ 1 : 3 ], beta$flair[ 1 : 3 ], 
                      lambda2$flair[ 1 : 3 ],
                      m, nu2$flair, 40L )
    flair_res <- flair_seg$parm[ c( 2 : 9 ), ]
    ## update delta
    
    x <- c( 0, 7, 10, 13, 13.5, 14 )
    y <- c( 0, 11, 15, 19, 17,  21 )
    spl <- smooth.spline( y, x )
    
    if( is.na( delta$flair[ 3 ] ) ) {
      shift_flair1 <- updateDelta3Flair( t1ce_image, 
                                         flair_data, flair_seg )
      shift_flair2 <- updateDelta12( flair_res[ 1, 3 ], 
                                     flair_res[ 2, 3 ],
                                     flair_res[ 1, 1 ], 
                                     flair_res[ 2, 1 ], 
                                     shift_flair1 )
      delta$flair[ c( 1, 2 ) ] <- delta$flair[ c( 1, 2 ) ] - 
        shift_flair2
      
      
      shift_flair1_ <- predict( spl, shift_flair1 )$y
      
      delta$flair[ 3 ] <- delta$flair[ 2 ] + 
        ( shift_flair1_ / delta_factor$flair ) ^ 2 / 2
      flair_shift <- c( delta$flair[ 1 : 2 ], shift_flair1, 
                        delta$flair[ 4 ] )
    } else {
      flair_shift <- delta$flair
      delta$flair[ 3 ] <- delta$flair[ 2 ] + 
        ( predict( spl, flair_shift[ 3 ] )$y /
            delta_factor$flair ) ^ 2 / 2
    }
    ## Export normalized images
    flair_intst <- flair_data$flair
    out_flair_norm <- gsub( ".nii.gz", "", out_flair_norm )
    writeNIfTI( nifti( flair_intst, datatype = 64 ) ,
                out_flair_norm,
                gzipped = T )
    save( flair_res, flair_shift, file = out_flair_data )
    ## update beta
    sigma2 <- flair_seg$parm[ 3, ]
    beta$flair[ 1 : 3 ] <- ( alpha$flair[ 1 : 3 ] + 1 ) * ( sigma2 )
    flair_data <- split4( flair_data$flair, flair_seg, flair_shift[ 3 ] )
    ## update m
    m <- flair_data$m
    b <- getB( m, a )
    flair_model <- initEst( flair_data$label, flair_data$intst )
    flair_seg <- pred( flair_model, delta$flair,
                       gamma$flair, alpha$flair, beta$flair,
                       lambda2$flair, a, b, m, nu2$flair, 
                       maxit$flair )
    ## flair
    flair_image <- flair_seg$image
    flair_image[ is.na( flair_image ) ] <- 0L
    flair_image[ flair_image >= 1 ] <- 4L
    flair_image[ flair_image <= -4 ] <- 4L
    flair_image[ flair_image == -1 ] <- 1L
    flair_image[ flair_image == -2 ] <- 2L
    flair_image[ flair_image == -3 ] <- 3L
    ## Get the initial results
    out_flair_seg <- gsub( ".nii.gz", "", out_flair_seg )
    writeNIfTI( nifti( flair_image, datatype = 2 ),
                filename = out_flair_seg, gzipped = TRUE )
    
  }
  
  out_t2_seg <- gsub( "_flair.nii.gz", "_t2_seg.nii.gz", outfile )
  out_t2_norm <- gsub( "_flair.nii.gz", "_t2_norm.nii.gz", outfile )
  out_t2_data <- gsub( "_flair.nii.gz", "_t2.RData", outfile )
  if( ( ! ( file.exists( out_t2_seg ) & 
          file.exists( out_t2_norm ) &
          file.exists( out_t2_data ) ) ) || 
      redo ) {
    ## split t2 to grey matter and white matter
    prop_bright <- propBright( t1ce_image, flair_image, images$t2 )
    t2_data <- splitT22( prop_bright, images$t2, t1ce_image, flair_image )
    m <- t2_data$m
    t2_model <- initEst( t2_data$label, t2_data$t2 )
    ## estimate parameters of t2 images without tumor 
    ## and CSF & necrosis
    # sink( '/home/hzhang/Documents/t2_output.txt' )
    t2_seg <- est( t2_model, delta$t2, gamma$t2, 
                   alpha$t2[ 1 : 2 ], beta$t2[ 1 : 2 ], 
                   lambda2$t2[ 1 : 2 ],
                   m, nu2$t2, 40L )
    t2_res <- t2_seg$parm[ c( 2 : 9 ), ]
    # sink()
    ## update delta
    if( is.na( delta$t2[ 3 ] ) ) {
      shift_t21 <- updateDelta3T2( prop_bright,
                                   t1ce_image, flair_image,
                                   t2_data, t2_seg ) 
      shift_t22 <- updateDelta12( t2_res[ 1, 2 ],
                                  t2_res[ 2, 2 ],
                                  t2_res[ 1, 1 ],
                                  t2_res[ 2, 1 ],
                                  shift_t21 )
      delta$t2[ 1 ] <- delta$t2[ 1 ] - shift_t22
      delta$t2[ 3 ] <- delta$t2[ 1 ] + 
        ( shift_t21 / delta_factor$t2 ) ^ 2 / 2
      t2_shift <- c( delta$t2[ 1 : 2 ], shift_t21, 
                     delta$t2[ 4 ] )
    } else {
      t2_shift <- delta$t2
      delta$t2[ 3 ] <- delta$t2[ 1 ] + 
        ( t2_shift[ 3 ] / delta_factor$t2 ) ^ 2 / 2
    }
    ## Export normalized images
    t2_intst <- t2_data$t2
    out_t2_norm <- gsub( ".nii.gz", "", out_t2_norm )
    save( t2_res, t2_shift, file = out_t2_data )
    ## update beta
    sigma2 <- t2_seg$parm[ 3, ]
    beta$t2[ 1 : 2 ] <- ( alpha$t2[ 1 : 2 ] + 1 ) * sigma2
    t2_data <- split3( t1ce_image, flair_image, prop_bright,
                       t2_data$t2, t2_seg, t2_shift[ 3 ] )
    ## update m
    m <- t2_data$m
    b <- getB( m, a )
    t2_model <- initEst( t2_data$label, t2_data$intst )
    t2_seg <- pred( t2_model, delta$t2, gamma$t2, alpha$t2, beta$t2,
                    lambda2$t2, a, b, m, nu2$t2, maxit$t2 )
    t2_seg$image[ t2_data$csf ] <- -4L
    ## t2
    t2_image <- t2_seg$image
    t2_image[ is.na( t2_image ) ] <- 0L
    t2_image[ t2_image >= 1 ] <- 4L
    t2_image[ t2_image <= -4 ] <- 4L
    t2_image[ t2_image == -1 ] <- 1L
    t2_image[ t2_image == -2 ] <- 2L
    ## Get the initial results
    out_t2_seg <- gsub( ".nii.gz", "", out_t2_seg )
    writeNIfTI( nifti( t2_image, datatype = 2 ),
                filename = out_t2_seg, gzipped = TRUE )
    writeNIfTI( nifti( t2_intst, datatype = 64 ) ,
                out_t2_norm,
                gzipped = T )
  }
}