#-------------------------------------------------------------------------------
#
#          MCMC
#
#-------------------------------------------------------------------------------

#' MCMC sampling
#'
#' MCMC sampling
#'
#' @param burnin the number of MCMC for burnin
#' @param niter the number of MCMC iterations
#' @param thin thinning of the MCMC chain
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
#' @param nu coefficient of variation for lognormal model
#' @param Ni yearly time series random effects
#' @param WdiPup daily time series random effects during pupping period
#' @param WdiMolt daily time series random effects during molting period
#' @param Xpg dataset for pup ground counts during pupping period 
#' @param Ypg dataset for nonpup ground counts during pupping period
#' @param Xpa dataset for pup aerial survey estimates during pupping period 
#' @param Ypa dataset for nonpup aerial survey estimates during pupping period
#' @param Ymg dataset for nonpup ground counts during molting period
#' @param Yma dataset for nonpup aerial survey estimates during molting period
#'
#' @return a list of the MCMC chains values
#'
#' @author Jay Ver Hoef
#' @export

MCMC <- function(burnin = 1, niter = 1000, thin = 2, eta, Ni,
  tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup,
  phiDyMolt, sigmaDyMolt, nu, WdiPup, WdiMolt, 
  Xpg, Ypg, Xpa, Ypa, Ymg, Yma)
{

  pdayrMin = min(Xpg$yday,Ypg$yday,Xpa$yday,Ypa$yday)
  pdayrMax = max(Xpg$yday,Ypg$yday,Xpa$yday,Ypa$yday)
  mdayrMin = min(Ymg$yday,Yma$yday)
  mdayrMax = max(Ymg$yday,Yma$yday)

loop1 = function(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
{
    for(i in 2:length(Ni)) {
      Ni.try = Ni
      U <- log(runif(1))
      Ni.try[i] <- rnorm(1, Ni[i], .1)
      LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni.try, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
        LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
      if(LLdif > U) Ni[i] <- Ni.try[i]
    }
   Ni
}
 
loop2 = function(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
{
		nyr = length(WdiPup[1, ])
	  ndayPup = length(WdiPup[, 1])
    for(j in 1:nyr){
      for(i in 2:ndayPup){
        U = log(runif(1))
        WdiPup.try = WdiPup
        WdiPup.try[i,j]  <- rnorm(1, WdiPup[i,j], .1)
        sumWdiPup.try = 0
        for(k in 1:nyr)
          sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup.try[,k], phiDyPup, sigmaDyPup)
        LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup.try, WdiMolt, sumWdiPup.try, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
          LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
        if(LLdif > U) {
          WdiPup <- WdiPup.try
          sumWdiPup <- sumWdiPup.try
        }
      }
    }
    list(WdiPup, sumWdiPup)
}

loop3 = function(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
{
		nyr = length(WdiMolt[1, ])
	  ndayMolt = length(WdiMolt[, 1])
		for(j in 1:nyr){
      for(i in 2:ndayMolt){
        WdiMolt.try = WdiMolt
        WdiMolt.try[i,j]  <- rnorm(1, WdiMolt[i,j], .1)
        sumWdiMolt.try = 0
        for(k in 1:nyr)
          sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt.try[,k], 
                           phiDyMolt, sigmaDyMolt)
        LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt.try, sumWdiPup, 
            sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
          LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
        if(LLdif > U) {
          WdiMolt <- WdiMolt.try
          sumWdiMolt <- sumWdiMolt.try
        }
      }
    }
    list(WdiMolt, sumWdiMolt)
}


	nyr = length(WdiPup[1, ])
	ndayPup = length(WdiPup[, 1])
	ndayMolt = length(WdiMolt[, 1])
  sumWdiPup = 0
  for(i in 1:nyr)
    sumWdiPup = sumWdiPup + dmvnormAR1(WdiPup[,i], phiDyPup, sigmaDyPup)
  sumWdiMolt = 0
  for(i in 1:nyr)
    sumWdiMolt = sumWdiMolt + dmvnormAR1(WdiMolt[,i], phiDyMolt, sigmaDyMolt)

  # ----------------------------------------------
  #                BURNIN
  # ----------------------------------------------

    if(burnin > 0) {
    for(iter in 1:burnin) {
    U <- log(runif(1))
    eta.try <- rnorm(1, eta, .02)
    LLdif <- LL(eta.try, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) eta <- eta.try
    U <- log(runif(1))
    tau.try <- rnorm(1, tau, .02)
    LLdif <- LL(eta, tau.try, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) tau <- tau.try
    U <- log(runif(1))
    xi.try <- rnorm(1, xi, .02)
    LLdif <- LL(eta, tau, xi.try, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) xi <- xi.try
    U <- log(runif(1))
    tauxi.try <- rnorm(1, tauxi, .02)
    LLdif <- LL(eta, tau, xi, tauxi.try, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) tauxi <- tauxi.try
    U <- log(runif(1))
    ppi.try <- rnorm(1, ppi, .02)
    LLdif <- LL(eta, tau, xi, tauxi, ppi.try, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) ppi <- ppi.try
    U <- log(runif(1))
    phiYr.try <- min(max(phiYr + runif(1)*0.04 - 0.02, 0.00001),.99999)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr.try, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) phiYr <- phiYr.try
    U <- log(runif(1))
    sigmaYr.try <- min(max(sigmaYr + runif(1)*.2 - .1, 0.01),10)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr.try, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) sigmaYr <- sigmaYr.try
    U <- log(runif(1))
    phiDyPup.try <- min(max(phiDyPup + runif(1)*0.04 - 0.02, 0.00001),.99999)
    sumWdiPup.try = 0
    for(i in 1:nyr)
      sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup[,i], phiDyPup.try, sigmaDyPup)    
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup.try, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup.try, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      phiDyPup <- phiDyPup.try
      sumWdiPup = sumWdiPup.try
    }
    U <- log(runif(1))
    sigmaDyPup.try <- min(max(sigmaDyPup + runif(1)*.2 - .1, 0.01),10)
    sumWdiPup.try = 0
    for(i in 1:nyr)
      sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup[,i], phiDyPup, sigmaDyPup.try)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup.try, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup.try, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      sigmaDyPup <- sigmaDyPup.try
      sumWdiPup = sumWdiPup.try
    }
    U <- log(runif(1))
    phiDyMolt.try <- min(max(phiDyMolt + runif(1)*0.04 - 0.02, 0.00001),.99999)
    sumWdiMolt.try = 0
    for(i in 1:nyr)
      sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt[,i], phiDyMolt.try, sigmaDyMolt)    
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt.try, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      phiDyMolt <- phiDyMolt.try
      sumWdiMolt = sumWdiMolt.try
    }
    U <- log(runif(1))
    sigmaDyMolt.try <- min(max(sigmaDyMolt + runif(1)*.2 - .1, 0.01),10)
    sumWdiMolt.try = 0
    for(i in 1:nyr)
      sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt[,i], phiDyMolt, sigmaDyMolt.try)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt.try, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      sigmaDyMolt <- sigmaDyMolt.try
      sumWdiMolt = sumWdiMolt.try
    }
    U <- log(runif(1))
    nu.try <- min(max(nu + runif(1)*.2 - .1, 0.01),10)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu.try, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) nu <- nu.try
    for(i in 2:length(Ni)) {
      Ni.try = Ni
      U <- log(runif(1))
      Ni.try[i] <- rnorm(1, Ni[i], .1)
      LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni.try, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
        LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
      if(LLdif > U) Ni[i] <- Ni.try[i]
    }
    i = 1
    j = 1
    for(j in 1:nyr){
      for(i in 2:ndayPup){
        U = log(runif(1))
        WdiPup.try = WdiPup
        WdiPup.try[i,j]  <- rnorm(1, WdiPup[i,j], .1)
        sumWdiPup.try = 0
        for(k in 1:nyr)
          sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup.try[,k], phiDyPup, sigmaDyPup)
        LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup.try, WdiMolt, sumWdiPup.try, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
          LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
        if(LLdif > U) {
          WdiPup <- WdiPup.try
          sumWdiPup <- sumWdiPup.try
        }
      }
    }
    for(j in 1:nyr){
      for(i in 2:ndayMolt){
        WdiMolt.try = WdiMolt
        WdiMolt.try[i,j]  <- rnorm(1, WdiMolt[i,j], .1)
        sumWdiMolt.try = 0
        for(k in 1:nyr)
          sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt.try[,k], 
                           phiDyMolt, sigmaDyMolt)
        LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt.try, sumWdiPup, 
            sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
          LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
            phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
            sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
        if(LLdif > U) {
          WdiMolt <- WdiMolt.try
          sumWdiMolt <- sumWdiMolt.try
        }
      }
    }
    cat("\r", "Burnin Number: ", iter)
  }
  cat("\n")
  }
  # ----------------------------------------------
  #          ITERATIONS WITH THINNING
  # ----------------------------------------------

  M <- vector("list", 16)
  M[[15]]=list()
  M[[16]]=list()
  iList = 1
  for(iter in 1:niter) {
    U <- log(runif(1))
    eta.try <- rnorm(1, eta, .02)
    LLdif <- LL(eta.try, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) eta <- eta.try
    U <- log(runif(1))
    tau.try <- rnorm(1, tau, .02)
    LLdif <- LL(eta, tau.try, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) tau <- tau.try
    U <- log(runif(1))
    xi.try <- rnorm(1, xi, .02)
    LLdif <- LL(eta, tau, xi.try, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) xi <- xi.try
    U <- log(runif(1))
    tauxi.try <- rnorm(1, tauxi, .02)
    LLdif <- LL(eta, tau, xi, tauxi.try, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) tauxi <- tauxi.try
    U <- log(runif(1))
    ppi.try <- rnorm(1, ppi, .02)
    LLdif <- LL(eta, tau, xi, tauxi, ppi.try, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) ppi <- ppi.try
    U <- log(runif(1))
    phiYr.try <- min(max(phiYr + runif(1)*0.04 - 0.02, 0.00001),.99999)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr.try, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) phiYr <- phiYr.try
    U <- log(runif(1))
    sigmaYr.try <- min(max(sigmaYr + runif(1)*.2 - .1, 0.01),10)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr.try, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) sigmaYr <- sigmaYr.try
    U <- log(runif(1))
    phiDyPup.try <- min(max(phiDyPup + runif(1)*0.04 - 0.02, 0.00001),.99999)
    sumWdiPup.try = 0
    for(i in 1:nyr)
      sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup[,i], phiDyPup.try, sigmaDyPup)	
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup.try, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup.try, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      phiDyPup <- phiDyPup.try
      sumWdiPup = sumWdiPup.try
    }
    U <- log(runif(1))
    sigmaDyPup.try <- min(max(sigmaDyPup + runif(1)*.2 - .1, 0.01),10)
    sumWdiPup.try = 0
    for(i in 1:nyr)
      sumWdiPup.try = sumWdiPup.try + dmvnormAR1(WdiPup[,i], phiDyPup, sigmaDyPup.try)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup.try, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup.try, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      sigmaDyPup <- sigmaDyPup.try
      sumWdiPup = sumWdiPup.try
    }
    U <- log(runif(1))
    phiDyMolt.try <- min(max(phiDyMolt + runif(1)*0.04 - 0.02, 0.00001),.99999)
    sumWdiMolt.try = 0
    for(i in 1:nyr)
      sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt[,i], phiDyMolt.try, sigmaDyMolt)    
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt.try, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      phiDyMolt <- phiDyMolt.try
      sumWdiMolt = sumWdiMolt.try
    }
    U <- log(runif(1))
    sigmaDyMolt.try <- min(max(sigmaDyMolt + runif(1)*.2 - .1, 0.01),10)
    sumWdiMolt.try = 0
    for(i in 1:nyr)
      sumWdiMolt.try = sumWdiMolt.try + dmvnormAR1(WdiMolt[,i], phiDyMolt, sigmaDyMolt.try)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt.try, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt.try, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) {
      sigmaDyMolt <- sigmaDyMolt.try
      sumWdiMolt = sumWdiMolt.try
    }
    U <- log(runif(1))
    nu.try <- min(max(nu + runif(1)*.2 - .1, 0.01),10)
    LLdif <- LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu.try, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin) -
      LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
        phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
        sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(LLdif > U) nu <- nu.try
    Ni = loop1(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    l2out = loop2(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
		WdiPup = l2out[[1]]
		sumWdiPup = l2out[[2]]
		l3out = loop3(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
		WdiMolt = l3out[[1]]
		sumWdiMolt = l3out[[2]]
    logLike = LL(eta, tau, xi, tauxi, ppi, phiYr, sigmaYr, phiDyPup, sigmaDyPup, 
          phiDyMolt, sigmaDyMolt, nu, Ni, WdiPup, WdiMolt, sumWdiPup, 
          sumWdiMolt, Xpg, Ypg, Xpa, Ypa, Ymg, Yma, pdayrMin, mdayrMin)
    if(iter%%thin == 0) {
      cat("\r", "Iter Number: ", iter)
      M[[1]] <- c(M[[1]],eta)
      M[[2]] <- c(M[[2]],tau)
      M[[3]] <- c(M[[3]],xi)
      M[[4]] <- c(M[[4]],tauxi)
      M[[5]] <- c(M[[5]],ppi)
      M[[6]] <- c(M[[6]],phiYr)
      M[[7]] <- c(M[[7]],sigmaYr)
      M[[8]] <- c(M[[8]],phiDyPup)
      M[[9]] <- c(M[[9]],sigmaDyPup)
      M[[10]] <- c(M[[10]],phiDyMolt)
      M[[11]] <- c(M[[11]],sigmaDyMolt)
      M[[12]] <- c(M[[12]],nu)
      M[[13]] <- cbind(M[[13]],Ni)
      M[[14]] <- c(M[[14]],logLike)
      M[[15]][[iList]] <- WdiPup
      M[[16]][[iList]] <- WdiMolt
      iList = iList + 1
    }
  }
  cat("\n")
  names(M) <- c("eta", "tau", "xi", "tauxi", "ppi", "phiYr", 
    "sigmaYr", "phiDyPup", "sigmaDyPup", "phiDyMolt", "sigmaDyMolt",
    "nu", "Ni", "logLike","WdiPup","WdiMolt")
  M
}

