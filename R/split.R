## split T1ce to CSF & necrosis, grey matter and white matter
splitT1ce3 <- function( t1ce, flair ) {
  
  label <- array( -4L, dim = dim( flair ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  q_flair <- quantile( flair, probs = .15, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .15, na.rm = T )
  
  ## find CSF & necrosis
  label[ flair < q_flair & t1ce < q_t1ce ] <- -1L
  ## Remove bright regions
  q_t1ce <- quantile( t1ce, probs = .90, na.rm = T )
  bright <- t1ce > q_t1ce
  label[ bright ] <- -4L
  
  tbd <- label == 0
  
  sub_flair <- flair[ tbd ]
  sub_t1ce <- t1ce[ tbd ]
  
  q_flair <- quantile( sub_flair, probs = c( .40, .50 ), 
                       na.rm =  T )
  q_t1ce <- quantile( sub_t1ce, probs =  c( .50, .60 ),
                      na.rm =  T )
  ## find white matter
  label[ tbd ][ sub_flair < q_flair[ 1 ] & sub_t1ce > q_t1ce[ 2 ] ] <- -3L
  ## find grey matter
  label[ tbd ][ sub_flair > q_flair[ 2 ] & sub_t1ce < q_t1ce[ 1 ] ] <- -2L
  
  label[ label == -4L ] <- NA_integer_
  
  ## normalize intensity
  mean_csf <- mean( t1ce[ label == -1 ], na.rm = T )
  mean_grey <- mean( t1ce[ label == -2 ], na.rm = T )
  mean_white <- mean( t1ce[ label == - 3 ], na.rm = T )
  
  t1ce <- ( t1ce - mean_csf ) / ( mean_white - mean_csf ) 
  mean_grey <- ( mean_grey - mean_csf ) / ( mean_white - mean_csf ) 
  res <- list( label = label, t1ce = t1ce,
               m = c( 0, mean_grey, 1 ) )
  return( res )
}
## split flair to CSF & necrosis, grey matter and white matter
splitFlair3 <- function( flair, t1ce_image ) {
  
  label <- array( -4L, dim = dim( flair ) )
  label[ ! is.nan( flair ) ] <- 0L
  
  ## Find the brightest 30% percent
  q_flair <- quantile( flair, probs = .70, na.rm = T )
  bright <- flair > q_flair
  
  # ## Find threshold dynamically
  # tumor <- flair[ flair > q_flair &
  #                 ( t1ce_image == 4 |
  #                   t1ce_image == 2 ) ] 
  # tumor <- tumor[ ! is.na( tumor ) ]
  # start_a <- quantile( tumor, probs = seq( 0, 1, length.out = 9 ) )
  # x_a <- kmeans( tumor, start_a )$centers
  # y_a <- order( x_a )
  # m_a <- x_a[ y_a ][ 3 ]
  # start_b <- quantile( tumor, probs = seq( 0, 1, length.out = 7 ) )
  # x_b <- kmeans( tumor, start_b )$centers
  # y_b <- order( x_b )
  # m_b <- x_b[ y_b ][ 3 ]
  # m_ab <- c( m_a, m_b )
  # m_ab <- m_ab[ order( m_ab ) ]
  # # m_m <- median( tumor[ tumor < m_ab[ 2 ] & tumor > m_ab[ 1 ] ], 
  # #              na.rm = T )
  # # m <- m_ab[ 2 ] - ( m_m - m_ab[ 1 ] )
  # m <- sum( m_ab ) / 2
  # bright <- flair > m
  label[ bright ] <- -4L
  
  ## Find CSF & necrosis
  q_flair <- quantile( flair, probs = .30, na.rm = T )
  cn <- t1ce_image == 1 & flair < q_flair
  label[ cn ] <- -1L
  tbd <- label == 0
  sub_flair <- flair[ tbd ]
  sub_image <- t1ce_image[ tbd ]
  
  ## Find white matter
  q_flair <- quantile( sub_flair, probs =  c( .40, .50 ), na.rm = T )
  wm <- sub_image == 3L & sub_flair < q_flair[ 1 ]
  label[ tbd ][ wm ] <- -2L
  
  ## Find grey matter
  gm <- sub_image == 2L & sub_flair > q_flair[ 2 ]
  label[ tbd ][ gm ] <- -3L
  
  label[ label == -4L ] <- NA_integer_
  
  ## Normalize intensity
  mean_csf <- mean( flair[ label == -1 ], na.rm =  T )
  mean_white <- mean( flair[ label == -2 ], na.rm = T )
  mean_grey <- mean( flair[ label == -3 ], na.rm = T )
  
  flair <- ( flair - mean_csf ) / ( mean_grey - mean_csf )
  mean_white <- ( mean_white - mean_csf ) / ( mean_grey - mean_csf )
  res <- list( label = label, flair = flair, 
               m = c( 0, mean_white, 1 ) )
  return( res )
}
## split t2 to grey matter and white matter
splitT22 <- function( prop_bright, t2, t1ce_image, flair_image ) {
  label <- array( -4L, dim = dim( t2 ) )
  label[ ! is.nan( t2 ) ] <- 0L
  
  ## Find CSF & necrosis and Tumor (brightest prop_bright)
  q_t2 <- quantile( t2, probs = 1 - prop_bright, na.rm = T )
  bright_t2 <- t2[ t2 > q_t2 & ( t1ce_image == 1 | 
                                 t1ce_image == 4 |
                                 flair_image == 1 |
                                 flair_image == 4 ) ]
  bright_t2 <- bright_t2[ ! is.na( bright_t2 ) ]
  b_km <- min( bright_t2 )
  bright <- t2 > b_km |
    t1ce_image == 1 | 
    t1ce_image == 4 |
    flair_image == 1 |
    flair_image == 4
  
  label[ bright ] <- -4L
  tbd <- label == 0
  sub_t2 <- t2[ tbd ]
  sub_t1ce <- t1ce_image[ tbd ]
  
  ## Find white matter
  q_t2 <- quantile( sub_t2, probs = c( .4, .7 ), na.rm = T )
  wm <- sub_t2 < q_t2[ 1 ]
  label[ tbd ][ wm ] <- -1L
  
  ## Find grey matter 
  gm <- sub_t1ce == 2L & sub_t2 > q_t2[ 2 ]
  label[ tbd ][ gm ] <- -2L
  
  label[ label == -4 ] <- NA_integer_
  
  ## normalize intensity
  mean_white <- mean( t2[ label == -1 ], na.rm = T )
  mean_grey <- mean( t2[ label == -2 ], na.rm = T )
  t2 <- ( t2 - mean_white ) / ( mean_grey - mean_white )
  res <- list( label = label, t2 = t2, m = c( 0, 1 ) )
  return( res )
}
## split t1ce or flair images for prediction
split4 <- function( x, x_seg, x_factor ) {
  label <- array( -4L, dim = dim( x ) )
  label[ ! is.nan( x ) ] <- 0L
  label[ x_seg$image == -1 ] <- -1L
  label[ x_seg$image == -2 ] <- -2L
  label[ x_seg$image == -3 ] <- -3L
  m_3 <- x_seg$parm[ 2, 3 ]
  sigma2_3 <- x_seg$parm[ 3, 3 ]
  ## threshold for enhancing tumor core
  m_4 <- m_3 + sqrt( sigma2_3 ) * x_factor
  label[ x > m_4 ] <- 4L
  m_1 <- x_seg$parm[ 2, 1 ]
  m_2 <- x_seg$parm[ 2, 2 ]
  
  label[ label == -4L ] <- NA_integer_
  
  res <- list( label = label, intst = x, 
               m = c( m_1, m_2, m_3, m_4 ) )
  return( res )
}
## split t2 images for prediction
split3 <- function( t1ce_image, flair_image, prop_bright,
                    x, x_seg, x_factor ) {
  label <- array( -4L, dim = dim( x ) )
  na_intst <- is.nan( x )
  label[ ! na_intst ] <- 0L
  t1ce_image[ na_intst ] <- 0L
  flair_image[ na_intst ] <- 0L
  
  q_t2 <- quantile( x, 1 - prop_bright, na.rm = T )
  x_sub <- x[ x > q_t2 ] 
  x_sub <- x_sub[ ! is.na( x_sub ) ]
  start <- quantile( x_sub, seq( 0, 1, length.out = 5 ) )
  mu <- kmeans( x_sub, start )$centers
  # csf <- ( t1ce_image == 1 |
  #          flair_image == 1 ) &
  #          x > mu[ 3 ]
  # csf <- which( csf )
  csf <- NULL
  
  label[ x_seg$image == -1 ] <- -1L
  label[ x_seg$image == -2 ] <- -2L
  m_2 <- x_seg$parm[ 2, 2 ]
  sigma2_2 <- x_seg$parm[ 3, 2 ]
  ## threshold for enhancing tumor core
  m_3 <- m_2 + sqrt( sigma2_2 ) * x_factor
  label[ x > m_3 ] <- 4L
  m_1 <- x_seg$parm[ 2, 1 ]
  
  label[ csf ] <- NA_integer_
  label[ label == -4L ] <- NA_integer_
  
  res <- list( label = label, intst = x, 
               m = c( m_1, m_2, m_3 ),
               csf = csf )
  return( res )
}
## Split edema || enh into to parts
## image: sub-image for each seperate tumor regions
splitFthrE <- function( image, t2_intst ) {
  post_mix_idx <- which( image == 2 | image == 4 )
  t2_mix <- t2_intst[ post_mix_idx ]
  valid_idx <- which( ! is.na( t2_mix ) )
  t2_mix <- t2_mix[ valid_idx ]
  post_mix_idx <- post_mix_idx[ valid_idx ]
  start <- c( min( t2_mix ), 
              max( t2_mix ) )
  clst <- kmeans( x = t2_mix, centers = start )$cluster
  label <- array( NA_integer_, dim = dim( t2_intst ) )
  mean1 <- mean( t2_mix[ clst == 1 ] )
  mean2 <- mean( t2_mix[ clst == 2 ] )
  if( mean1 > mean2 ) {
    clst[ clst == 1 ] <- -2L
    clst[ clst == 2 ] <- -1L
    m <- c( mean2, mean1 )
  } else {
    clst[ clst == 1 ] <- -1L
    clst[ clst == 2 ] <- -2L
    m <- c( mean1, mean2 )
  }
  label[ post_mix_idx ] <- clst
  res <- list( label = label, intst = t2_intst, 
               m = m )
  return( res )
}
## Split CSF into parts
splitFthrC <- function( image, t1ce_intst ) {
  post_mix_idx <- which( image == 5 | 
                         image == 6 |
                         image == 1 )
  t1ce_mix <- t1ce_intst[ post_mix_idx ]
  valid_idx <- which( ! is.na( t1ce_mix ) )
  t1ce_mix <- t1ce_mix[ valid_idx ]
  post_mix_idx <- post_mix_idx[ valid_idx ]
  
  n_clst <- 5
  start <- quantile( t1ce_mix, seq( 0, 1, length.out = n_clst ) )
  k_m <- kmeans( t1ce_mix, start )
  k_centers <- k_m$centers
  k_cluster <- k_m$cluster
  ord <- order( k_centers )
  valid_idx2 <- which( k_cluster != ord[ 1 ] )
  clst <- k_cluster[ valid_idx2 ]
  t1ce_mix <- t1ce_mix[ valid_idx2 ]
  post_mix_idx <- post_mix_idx[ valid_idx2 ]
  
  mean1 <- k_centers[ ord[ 2 ] ]
  mean2 <- k_centers[ ord[ n_clst - 1 ] ]
  m <- c( mean1, mean2 )
  clst[ clst == ord[ 2 ] ] <- -1L
  clst[ clst == ord[ 4 ] ] <- -2L
  clst[ clst != -1L & clst != -2L ] <- -0L
  label <- array( NA_integer_, dim = dim( t1ce_intst ) )
  label[ post_mix_idx ] <- clst
  res <- list( label = label, intst = t1ce_intst, 
               m = m )
  return( res )
}