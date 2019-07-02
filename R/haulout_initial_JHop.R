#-------------------------------------------------------------------------------
#
#          haulout_initial
#
#-------------------------------------------------------------------------------

#' Initialize haulout MCMC chain
#'
#' Initialize haulout MCMC chain
#'
#' @param season either "spring" or "fall".
#'
#' @return a list of inputs for stock_haulout_estimation() function
#'
#' @author Jay Ver Hoef
#' @export

haulout_initial_JHop = function(season)
{

  if(tolower(substr(season,1,4))=='spr') {
    data(dHOspr)
    dHO = dHOspr
  }
  if(tolower(substr(season,1,4))=='fal') {
    data(dHOfal)
    dHO = dHOfal
  }

  
  # make sure speno's are factors
  dHO$speno = as.factor(as.character(dHO$speno))
  HOIDs = levels(dHO$speno)
  #order by datadatetime within speno, and remove duplicate times
  dHO = dHO[order(dHO$speno, dHO$date_time),]

  # some initial beta parameters
  mu = 0.4
  b = .11
  a = b*mu/(1 - mu)
  X <- model.matrix(y ~ dystd + I(dystd^2), data = dHO)
  Z <- model.matrix(y ~ -1 + speno, data = dHO)
  beta = rep(0, times = dim(X)[2])
  Zi = as.integer(dHO$speno)

  # use beta distribution to get some initial parameters
  LLho = function(theta,y,X)
  {
    # theta[1] is b
    # rest of theta is part of linear model
    Xb = X %*% theta[2:length(theta)]
    mu = exp(Xb)/(1 + exp(Xb))
    a = mu*theta[1]
    b = (1-mu)*theta[1]
      -sum(dbeta(y,a,b, log = TRUE))
  }
  optout1 = optim(c(.5,.5,.5,-.5), LLho, y = dHO$y, X = X)


  # initialize MCMC parameters
  phi = optout1$par[1]
  beta = optout1$par[2:length(optout1$par)]
  alpha = .1
  sigGam = 1
  sigmaARexp = 1
  gam = rep(0, times = max(Zi))
  tdif = list()
  eps = list()
  i = 1
  for(i in 1:length(levels(dHO$speno))) {
    hrs = dHO[dHO$speno == levels(dHO$speno)[i],'yrhr0']
    tdif[[i]] = hrs[2:length(hrs)] -  hrs[1:(length(hrs) - 1)]
  }
  for(i in 1:length(tdif)){
    eps[[i]] = (dHO[dHO$speno == levels(dHO$speno)[i],'y'] - 0.5)*2
    eps[[i]][1] = 0
  }


  phi_tune = 0.1
  beta_tune = rep(0.02, times = length(beta))
  alpha_tune = .2
  sgam_tune = 0.2
  sAR1_tune = 0.01
  gam_tune = rep(.05, times = length(gam))
  eps_tune = rep(.02, times = length(gam))
  
  list(dHO = dHO, phi = phi, beta = beta, alpha = alpha, sigGam = sigGam, 
    sigmaARexp = sigmaARexp, gam = gam, eps = eps, y = dHO$y, X = X, 
    Zi = Zi, tdif = tdif, phi_tune = phi_tune, beta_tune = beta_tune,
    alpha_tune = alpha_tune, sgam_tune = sgam_tune, 
    sAR1_tune = sAR1_tune, gam_tune = gam_tune, eps_tune = eps_tune)
}

