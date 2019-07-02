#-------------------------------------------------------------------------------
#
#          LL
#
#-------------------------------------------------------------------------------

#' loglikelihood for model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param eta overall intercept
#' @param tau effect for survey method being ground count
#' @param xi effect for surveys in the pupping period
#' @param tauxi ground count by pupping period interaction
#' @param ppi effect for pup counts
#' @param phiYr AR1 autocorrelation parameter for yearly time series
#' @param sigmaYr AR1 variance parameter for yearly time series
#' @param phiDyPup AR1 autocorrelation parameter for daily time series during pupping period
#' @param sigmaDyPup AR1 variance parameter for daily time series during pupping period
#' @param phiDyMolt AR1 autocorrelation parameter for daily time series during molting period
#' @param sigmaDyMolt  variance parameter for daily time series during molting period
#' @param sigmaZ variance parameter for independent residuals
#' @param Ni yearly time series random effects
#' @param WdiPup daily time series random effects during pupping period
#' @param WdiMolt daily time series random effects during molting period
#' @param sumWdiPup sum of the likelihood of daily random effects during pupping period
#' @param sumWdiMolt sum of the likelihood of daily random effects during the molting period
#' @param Xpg dataset for pup ground counts during pupping period 
#' @param Ypg dataset for nonpup ground counts during pupping period
#' @param Xpa dataset for pup aerial survey estimates during pupping period 
#' @param Ypa dataset for nonpup aerial survey estimates during pupping period
#' @param Ymg dataset for nonpup ground counts during molting period
#' @param Yma dataset for nonpup aerial survey estimates during molting period
#' @param pdayrMin minimum yday for all data sets during pupping period
#' @param mdayrMin minimum yday for all data sets during molting period
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef
#' @export

LL <- function(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup,
  phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, sumWdiMolt, 
  Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
{
  sigmaZ = sqrt(log(nu^2 + 1))
  logMuXpg = eta + Ni[Xpg$year - 1991] + tau + xi + tauxi + ppi +
    WdiPup[cbind(Xpg$yday - pdayrMin + 1, Xpg$year - 1991)] - sigmaZ^2/2
  logMuYpg = eta + Ni[Ypg$year - 1991] + tau + xi + tauxi +
    WdiPup[cbind(Ypg$yday - pdayrMin + 1, Ypg$year - 1991)] - sigmaZ^2/2
  logMuXpa = eta + Ni[Xpa$year - 1991] + xi + ppi +
    WdiPup[cbind(Xpa$yday - pdayrMin + 1, Xpa$year - 1991)] - sigmaZ^2/2
  logMuYpa = eta + Ni[Ypa$year - 1991] + xi +
    WdiPup[cbind(Ypa$yday - pdayrMin + 1, Ypa$year - 1991)] - sigmaZ^2/2
  logMuYmg = eta + Ni[Ymg$year - 1991] + tau +
    WdiMolt[cbind(Ymg$yday - mdayrMin + 1, Ymg$year - 1991)] - sigmaZ^2/2
  logMuYma = eta + Ni[Yma$year - 1991] +
    WdiMolt[cbind(Yma$yday - mdayrMin + 1, Yma$year - 1991)] - sigmaZ^2/2
  logSdXpg = sigmaZ
  logSdYpg = sigmaZ
  logSdXpa = sqrt(log((Xpa$pupSE/Xpa$pupEst)^2 + 1))
  logSdYpa = sqrt(log((Ypa$aduSE/Ypa$aduEst)^2 + 1))
  logSdYmg = sigmaZ
  logSdYma = sqrt(log((Yma$aduSE/Yma$aduEst)^2 + 1))
  sum(dlnorm(Xpg$pup, meanlog = logMuXpg, sdlog = logSdXpg, 
    log = TRUE)) + 
  sum(dlnorm(Ypg$adu, meanlog = logMuYpg, sdlog = logSdYpg, 
    log = TRUE)) +
  sum(dlnorm(Xpa$pupEst, meanlog = logMuXpa, sdlog = logSdXpa, 
    log = TRUE)) + 
  sum(dlnorm(Ypa$aduEst, meanlog = logMuYpa, sdlog = logSdYpa, 
    log = TRUE)) + 
  sum(dlnorm(Ymg$adu, meanlog = logMuYmg, sdlog = logSdYmg, 
    log = TRUE)) +
  sum(dlnorm(Yma$aduEst, meanlog = logMuYma, sdlog = logSdYma, 
    log = TRUE)) +
  dmvnormAR1(Ni, phiYr, sigmaYr) +
  sumWdiPup + sumWdiMolt + 
	  dnorm(sigmaDyPup,.5,.1) +
		dnorm(sigmaDyMolt,.5,.1)
}

