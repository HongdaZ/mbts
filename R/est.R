## estimate parameters of t1ce, flair images without tumor or
## t2 images without tumor and CSF & necrosis
est <- function( model, delta, gamma, 
                  alpha, beta, lambda2, m, nu2, maxit ) {
  .Call( "est", model, delta, gamma, 
         alpha, beta, lambda2, m, nu2, maxit )
}