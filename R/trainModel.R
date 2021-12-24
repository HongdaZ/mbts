# # Train the model
# # beta and nu2 influenced by normalization
# trainModel <- function( patient, delta = 5 ^ 2, gamma = 1, 
#                         alpha = rep( 10, 4 ),
#                         beta = rep( 1, 4 ),
#                         lambda2 = rep( 1 / ( 36 * 4 ) , 4 ), 
#                         a = 5,
#                         nu2 = rep( .25, 3 ), 
#                         maxit = 50L ) {
#   
#   l_intst <- splitCWG( patient ) # splitCWGX for testing
#   flair_model <- initTrn( l_intst$label_flair, l_intst$intst_flair, 
#                           "flair" )
#   # -1, -2, -3, > 0
#   m <- priorMode( flair_model )
#   b <- getB( m, a )
#   flair_res <- estParm( flair_model, delta, gamma,
#            alpha, beta, lambda2, a, b, m, nu2, 1L )
#   return( flair_res )
# }
