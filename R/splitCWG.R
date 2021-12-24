# # Split CSF & necrosis, white matter and grey matter in training data
# # csf = 1, wm = 2, g = 3, TBD = 0, o.w. = NA 
# splitCWG <- function( patient ) {
#   
#   lbl_ <- readNIfTI( patient[ 2 ], reorient = FALSE )@.Data
#   flair_ <- readNIfTI( patient[ 1 ], reorient = FALSE )@.Data
#   t2_ <- readNIfTI( patient[ 5 ], reorient = FALSE )@.Data
#   t1ce_ <- readNIfTI( patient[ 4 ], reorient = FALSE )@.Data
# 
#   flair_[ flair_ == 0 ] <- NaN
#   t2_[ t2_ == 0 ] <- NaN
#   t1ce_[ t1ce_ == 0 ] <- NaN
# 
#   # non_valid <- lbl_ != 0
#   # 
#   # flair[ non_valid ] <-NaN
#   # t2[ non_valid ] <- NaN
#   # t1ce[ non_valid ] <- NaN
#   flair <- changeA( flair_, lbl_ )
#   t2 <- changeA( t2_, lbl_ )
#   t1ce <- changeA( t1ce_, lbl_ )
#   
#   # Find csf & necrosis
#   q_flair <- quantile( flair, probs = c( .20, .40, .50, .01 ), na.rm = T, 
#                        names = F )
#   q_t2 <- quantile( t2, probs = c( .80, .01 ), na.rm = T, names = F )
#   
#   # csf <- flair < q_flair[ 1 ] & t2 > q_t2[ 1 ]
#   # lbl_[ csf ] <- -1L
#   lbl_ <- changeB( lbl_, flair, q_flair[ 1 ], t2, q_t2[ 1 ], -1L )
#   
#   # Find white matter
#   q_t1ce <- quantile( t1ce, probs = c( .60, .50, .01 ), na.rm = T, names = F )
#   
#   # wm <- flair < q_flair[ 2 ] & t1ce > q_t1ce[ 1 ]
#   # lbl_[ wm ] <- -2L
#   lbl_ <- changeB( lbl_, flair, q_flair[ 2 ], t1ce, q_t1ce[ 1 ], -2L )
#   
#   # Find grey matter
#   # gm <- t1ce < q_t1ce[ 2 ] & flair > q_flair[ 3 ]
#   # lbl_[ gm ] <- -3L
#   lbl_ <- changeB( lbl_, t1ce, q_t1ce[ 2 ], flair, q_flair[ 3 ], -3L )
#   
#   # Remove the darkest 1% of the voxels
#   # lbl_flair <- lbl_
#   # lbl_flair[ flair < q_flair[ 4 ] ] <- NA_integer_
#   lbl_flair <- changeC( lbl_, flair, q_flair[ 4 ] )
#   
#   # lbl_t1ce <- lbl_
#   # lbl_t1ce[ t1ce < q_t1ce[ 3 ] ] <- NA_integer_
#   lbl_t1ce <- changeC( lbl_, t1ce, q_t1ce[ 3 ] )
#   
#   # lbl_t2 <- lbl_
#   # lbl_t2[ t2 < q_t2[ 2 ] ] <- NA_integer_
#   lbl_t2 <- changeC( lbl_, t2, q_t2[ 2 ] )
#   
#   
#   res <- list( label_flair = lbl_flair, label_t1ce = lbl_t1ce, 
#                label_t2 = lbl_t2, 
#                intst_flair = flair_, intst_t1ce = t1ce_,
#                intst_t2 = t2_ )
#   return( res )
# }