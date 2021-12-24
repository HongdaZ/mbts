# .Call interface doesn't invoke copy-on-modify sematics
call_interface <- function( x ) {
  .Call( "call_interface", x )
}