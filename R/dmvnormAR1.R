#-------------------------------------------------------------------------------
#
#          dmvnormAR1
#
#-------------------------------------------------------------------------------

#' loglikelihood for AR1 model
#'
#' evaluates to the loglikelihood for an AR1 model conditional on the first value being 0
#'
#' @param y the data
#' @param phi the autocorrelation parameter
#' @param sigmaAR1 the variance parameter
#
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

dmvnormAR1 = function(y, phi, sigmaAR1) {
  n = length(y)
  -(n - 1)/2*log(2*pi) - (n - 1)*log(sigmaAR1) -
  sum((y[2:n] - phi*y[1:(n - 1)])^2)/(2*sigmaAR1^2)
}


