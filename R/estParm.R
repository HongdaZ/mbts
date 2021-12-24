# fit model to estimate mu1, mu2, mu3 and sigma's
estParm <- function( model, delta, gamma, 
                            alpha, beta, lambda2, a, b, m, nu2, maxit ) {
  .Call( "estParm", model, delta, gamma, 
         alpha, beta, lambda2, a, b, m, nu2, maxit )
}