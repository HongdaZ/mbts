# predict the labels for T1ce or T2 images
pred <- function( model, delta, gamma,
                  alpha, beta, lambda2, a, b, m, nu2, maxit ) {
  .Call( "pred", model, delta, gamma, 
         alpha, beta, lambda2, a, b, m, nu2, maxit )
}