## Furtherly segment edema in FLAIR
estF <- function( model, delta, gamma, 
                 alpha, beta, lambda2, m, nu2, maxit ) {
  .Call( "estF", model, delta, gamma, 
         alpha, beta, lambda2, m, nu2, maxit )
}