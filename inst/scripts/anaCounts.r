library(JHop)

  # create the data
data(dgrnd)
data(dair)
dpup = dair[dair$age == 'pup',]
dadu = dair[dair$age == 'adu',]
dpup$date1 = paste0(dpup$year,',',dpup$yday)
dadu$date1 = paste0(dadu$year,',',dadu$yday)
dair1 = merge(dadu,dpup, by = 'date1', all.x = TRUE)
dair1 = dair1[,c('year.x', 'moday.x', 'yday.x', 'PopEst.x',
  'PopSE.x', 'PopEst.y', 'PopSE.y')]
names(dair1) = c('year','moday','yday','aduEst','aduSE','pupEst','pupSE')
dair1$season = 'pup'
dair1[dair1$yday > 210,'season'] = 'molt'
dgrnd$season = 'pup'
dgrnd[dgrnd$yday > 210,'season'] = 'molt'
Xpg = dgrnd[dgrnd$season == 'pup' & !is.na(dgrnd$pup),]
Xpg[Xpg$pup == 0,'pup'] = 0.01
Ypg = dgrnd[dgrnd$season == 'pup',]
Xpa = dair1[dair1$season == 'pup',]
Ypa = dair1[dair1$season == 'pup',]
Ymg = dgrnd[dgrnd$season == 'molt',]
Yma = dair1[dair1$season == 'molt',]

  eta = 8.2
  Ni = c(0,  0.07, 0.11,  0.17,  0.29,
    0.16, -0.18, -0.17, -0.31, -0.37,
    -0.51, -0.26, -0.28, -0.35, -0.54,
    -0.79, -0.97, -0.88, -0.64, -0.59,
    -0.69, -0.42, -.5, -.5, -.5, 
    -.5)
  tau = -.3
  xi = .05
  tauxi = -.3
  ppi = -1.3
  phiYr = 0.8
  sigmaYr = 0.2
  phiDyPup = .94
  sigmaDyPup = 0.5
  phiDyMolt = .7
  sigmaDyMolt = 0.5
  nu = 0.46
  pdayrMin = min(Xpg$yday, Ypg$yday, Xpa$yday, Ypa$yday)
  pdayrMax = max(Xpg$yday, Ypg$yday, Xpa$yday, Ypa$yday)
  mdayrMin = min(Ymg$yday, Yma$yday)
  mdayrMax = max(Ymg$yday, Yma$yday)
  WdiPup = matrix(0, ncol = length(Ni), nrow = pdayrMax - pdayrMin + 1)
  WdiMolt = matrix(0, ncol = length(Ni), nrow = mdayrMax - mdayrMin + 1)

#undebug(MCMC)
#undebug(LL)
MCMCout = MCMC(burnin = 0, niter = 4000, thin = 4, eta = eta, Ni = Ni,
  tau = tau, xi = xi, tauxi = tauxi, ppi = ppi, phiYr = phiYr, 
  sigmaYr = sigmaYr, phiDyPup = phiDyPup, sigmaDyPup = sigmaDyPup,
  phiDyMolt = phiDyMolt, sigmaDyMolt = sigmaDyMolt, nu = nu,
  WdiPup = WdiPup, WdiMolt = WdiMolt, Xpg = Xpg, Ypg = Ypg, 
  Xpa = Xpa, Ypa = Ypa, Ymg = Ymg, Yma = Yma)
M0 = MCMCout
nM = length(M0[['eta']])
save(M0, file = paste0('/media/jay/Hitachi2GB/00NMML/ActiveRPack/',
  'Harb_Seal_Surv/JHop/data/M0.rda'))


M = M0
plot(M[['eta']], type = 'l')
plot(M[['tau']], type = 'l')
plot(M[['xi']], type = 'l')
plot(M[['tauxi']], type = 'l')
plot(M[['ppi']], type = 'l')
plot(M[['phiYr']], type = 'l')
plot(M[['sigmaYr']], type = 'l')
plot(M[['phiDyPup']], type = 'l')
plot(M[['sigmaDyPup']], type = 'l')
plot(M[['phiDyMolt']], type = 'l')
plot(M[['sigmaDyMolt']], type = 'l')
plot(M[['nu']], type = 'l')
plot(M[['Ni']][2,], type = 'l')
plot(M[['logLike']], type = 'l')
junk = M[['WdiPup']]
plot(junk[[1000]][,20], type = 'l')

# data(M)
# nM = length(M[['eta']])
# M0 = M
MCMCout1 = MCMC(burnin = 0, niter = 20000, thin = 20, 
  eta = M0[['eta']][nM], 
  Ni = M0[['Ni']][,nM], 
  tau = M0[['tau']][nM],
  xi = M0[['xi']][nM], 
  tauxi = M0[['tauxi']][nM], 
  ppi = M0[['ppi']][nM], 
  phiYr = M0[['phiYr']][nM], 
  sigmaYr = M0[['sigmaYr']][nM], 
  phiDyPup = M0[['phiDyPup']][nM], 
  sigmaDyPup = M0[['sigmaDyPup']][nM], 
  phiDyMolt = M0[['phiDyMolt']][nM], 
  sigmaDyMolt = M0[['sigmaDyMolt']][nM], 
  nu = M0[['nu']][nM],
  WdiPup = M0[['WdiPup']][[nM]], 
  WdiMolt = M0[['WdiMolt']][[nM]], 
  Xpg = Xpg, Ypg = Ypg, Xpa = Xpa, Ypa = Ypa, Ymg = Ymg, Yma = Yma)

plot(MCMCout1[['eta']], type = 'l')
plot(MCMCout1[['tau']], type = 'l')
plot(MCMCout1[['xi']], type = 'l')
plot(MCMCout1[['tauxi']], type = 'l')
plot(MCMCout1[['ppi']], type = 'l')
plot(MCMCout1[['phiYr']], type = 'l')
plot(MCMCout1[['sigmaYr']], type = 'l')
plot(MCMCout1[['phiDyPup']], type = 'l')
plot(MCMCout1[['sigmaDyPup']], type = 'l')
plot(MCMCout1[['phiDyMolt']], type = 'l')
plot(MCMCout1[['sigmaDyMolt']], type = 'l')
plot(MCMCout1[['nu']], type = 'l')
plot(MCMCout1[['Ni']][2,], type = 'l')

# M0 = MCMCout1
# nM = length(M0[['eta']])
M = MCMCout1
nM = length(M[['eta']])
save(M, file = paste0('/media/jay/Hitachi2GB/00NMML/ActiveRPack/',
  'Harb_Seal_Surv/JHop/data/M.rda'))


