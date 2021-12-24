# lbl [ flair < q_flair[ 1 ] & t2 > q_t2[ 1 ] ] <- k
changeB <- function( label, image1, q1, image2, q2, k ) {
  .Call( "changeB", label, image1, q1, image2, q2, k )
}