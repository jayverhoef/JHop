figpath = paste0('/media/jay/Hitachi2GB/00NMML/ActiveRPack/',
	'Harb_Seal_Surv/JHop/inst/figures/')

## ----loadLibrary, echo = FALSE, include = FALSE--------------------------
  # load the library anew each time
  library(JHop)
  library(sp)

## ----createData, echo = FALSE, include = FALSE, cache = TRUE-------------
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
data(M)


## ----makeHOdata, echo = FALSE, include = FALSE, cache = TRUE-------------
data(dHOfal)
data(dHOspr)
Xfal <- model.matrix(y ~ dystd + I(dystd^2), data = dHOfal)
Xspr <- model.matrix(y ~ dystd + I(dystd^2), data = dHOspr)
Zfal <- model.matrix(y ~ -1 + speno, data = dHOfal)
Zspr <- model.matrix(y ~ -1 + speno, data = dHOspr)
data(sHOe_GLAC)
MHOaug = sHOe_GLAC$M

## ----Tab-Abun, echo= FALSE, include = FALSE, cache = TRUE----------------
  library(xtable)
  NiMat = M[['Ni']]
  etaVec = M[['eta']]
  mnWdiMolt = 0
  for(i in 1:1000) mnWdiMolt = mnWdiMolt + M[['WdiMolt']][[i]][13,]
  mnWdiMolt = mnWdiMolt/1000
  beta0 = unlist(lapply(MHOaug[['beta']], function(x){x[[1]]}))
  AbunMat = matrix(0, nrow = 1000, ncol = length(NiMat[,1]))
  for(i in 1:1000)
    AbunMat[i,] = (exp(etaVec[i] + NiMat[,i] + mnWdiMolt))*
      (1+exp(beta0[i]))/exp(beta0[i])
    AbunTable = cbind(
      apply(AbunMat,2,quantile,probs=.025),
      apply(AbunMat,2,quantile,probs=.05), 
      apply(AbunMat,2,mean),
      apply(AbunMat,2,quantile,probs=.95),
      apply(AbunMat,2,quantile,probs=.975)
      )
  string = as.character(1992:2017)
  rownames(AbunTable) = string

## ----results = 'asis', echo = FALSE--------------------------------------
  print(
    xtable(AbunTable, 
      align = c('l',rep('l', times = length(AbunTable[1,]))),
      digits = c(0,0,0,0,0,0),
      caption = 'Abundance Table',
      label = 'tab:Abun'
    ),
    size = 'footnotesize',
    include.rownames = TRUE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  ) 

## ----Tab-trendProgBack, echo= FALSE, include = FALSE, cache = TRUE-------
  colnames(AbunMat) = 1992:2017
  linTrendMCMC = function(AbunMat,startyr,stopyr)
  {
    startindx = (1:dim(AbunMat)[2])[colnames(AbunMat) == as.character(startyr)]
    stopindx = (1:dim(AbunMat)[2])[colnames(AbunMat) == as.character(stopyr)]
    trendi = rep(0, times = dim(AbunMat)[1])
    for(i in 1:dim(AbunMat)[1]) {
      y = AbunMat[i,startindx:stopindx]
      x = 1:length(y)
      trendi[i] = coef(lm(y ~ x))[2]
    }
    trendi
  } 
  trendsProgBack = NULL
  for(j in 2013:1992) {
    trendi = linTrendMCMC(AbunMat, j, 2017)
    trendsProgBack = rbind(trendsProgBack,
      data.frame(years = paste0(j,'-2017'),
        low95 = quantile(trendi, probs = 0.025),
        low90 = quantile(trendi, probs = 0.05),
        Est = mean(trendi), 
        upp90 = quantile(trendi,probs = 0.95),
        upp95 = quantile(trendi,probs = 0.975)
      )
    )
  }

## ----results = 'asis', echo = FALSE--------------------------------------
    print(
    xtable(trendsProgBack, 
      align = c('l',rep('l', times = length(trendsProgBack[1,]))),
      digits = c(0,0,0,0,0,0,0),
      caption = 'Trends Progressively Backwards Table',
      label = 'tab:trendsProgBack'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  ) 

## --------------- TREND ESTIMATES, 10 YEAR INCREMENTS

  trends10YrInc = NULL
  for(j in 1992:2008) {
    trendi = linTrendMCMC(AbunMat, j, j+9)
    trends10YrInc = rbind(trends10YrInc,
      data.frame(years = paste0(j,'-',j+9),
        low95 = quantile(trendi, probs = 0.025),
        low90 = quantile(trendi, probs = 0.05),
        Est = mean(trendi), 
        upp90 = quantile(trendi,probs = 0.95),
        upp95 = quantile(trendi,probs = 0.975)
      )
    )
  }
  
## ----results = 'asis', echo = FALSE--------------------------------------
    print(
    xtable(trends10YrInc, 
      align = c('l',rep('l', times = length(trends10YrInc[1,]))),
      digits = c(0,0,0,0,0,0,0),
      caption = '10 Year Trends Incremented Yearly',
      label = 'tab:trends10YrInc'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  ) 
  
## ----grndAirScatters, echo=FALSE, include = FALSE, fig.height = 8, fig.width = 8, cache = TRUE----

png(filename = paste0(figpath,'grndAirScatters.png'), width = 550, height = 550) 
  dgrnd$date1 = paste0(dgrnd$year,',',dgrnd$yday)
  dair1$date1 = paste0(dair1$year,',',dair1$yday)
  dboth = merge(dgrnd,dair1, by = 'date1')
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  par(mar = c(5,5,5,1))
  plot(dboth[dboth$season.x == 'pup',c('adu','aduEst')], pch = 19,
    xlim = c(0,3000), ylim = c(0,3000),
    xlab = 'Ground Counts of Nonpups', cex.lab = 1.5, cex.axis = 1.5,
    ylab = 'Aerial Survey Estimates of Nonpups',
    main = 'Pupping Season', cex.main = 2)
  lines(c(0,3000),c(0,3000), lwd = 2)
  plot(dboth[dboth$season.x == 'pup',c('pup','pupEst')], pch = 19,
    xlim = c(0,1600), ylim = c(0,1600),
    xlab = 'Ground Counts of Pups', cex.lab = 1.5, cex.axis = 1.5,
    ylab = 'Aerial Survey Estimates of Pups',
    main = 'Pupping Season', cex.main = 2)
  lines(c(0,1600),c(0,1600), lwd = 2)
  plot(dboth[dboth$season.x == 'molt',c('adu','aduEst')], pch = 19,
    xlim = c(0,2100), ylim = c(0,2100),
    xlab = 'Ground Counts of Nonpups', cex.lab = 1.5, cex.axis = 1.5,
    ylab = 'Aerial Survey Estimates of Nonpups',
    main = 'Molting Season', cex.main = 2)
  lines(c(0,2100),c(0,2100), lwd = 2)
dev.off()

## ----postProp, echo=FALSE, include = FALSE, fig.height = 10, fig.width = 10, cache = TRUE----

png(filename = paste0(figpath,'postProp.png'), width = 680, height = 1020) 
  layout(matrix(1:6, ncol = 2, byrow = TRUE))
  par(mar = c(6,6,6,2))
  plot(density(exp(M[['xi']]), bw = .07),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Season Main Effect',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(a)', cex = 2.5, adj = -.25, padj = -1.4)
  plot(density(exp(M[['tau']]), bw = .04),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Survey Method Main Effect',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(b)', cex = 2.5, adj = -.25, padj = -1.4)
  plot(density(exp(M[['tauxi']]), bw = .04),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Season by Method Interaction',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(c)', cex = 2.5, adj = -.25, padj = -1.4)
  plot(density(exp(M[['ppi']]), bw = .002),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Pup Proportion',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(d)', cex = 2.5, adj = -.25, padj = -1.4)
  plot(density(exp(M[['xi']]) + exp(M[['xi']] + M[['ppi']]), bw = .1),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Total Aerial Estimates During Pupping',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(e)', cex = 2.5, adj = -.25, padj = -1.4)
  plot(density(exp(M[['tau']] + M[['xi']] + M[['tauxi']]) + 
		exp(M[['tau']] + M[['xi']] + M[['tauxi']] + M[['ppi']]), bw = .07),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Total Ground Counts During Pupping',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)
  mtext('(f)', cex = 2.5, adj = -.25, padj = -1.4)
dev.off()




plot(density(exp(M[['tau']] + M[['xi']] + M[['tauxi']]) + exp(M[['tau']] + M[['xi']] + M[['tauxi']] + M[['ppi']]), bw = .07),
    xlab = 'Proportional Effect',
    ylab = 'Posterior Density', lwd = 3,
    main = 'Aerial Counts During Pupping',
    cex.lab = 2, cex.axis = 1.5, cex.main = 1.7)

## ----trendFits, echo=FALSE, include = FALSE, fig.height = 9, fig.width = 15, cache = TRUE----

RGBcolRaw = c(189,189,189) 
RGBcolEst = c(254,217,118)
RGBcolPup = c(55,126,184) 
RGBcolAdu = c(77,175,74)
RGBcolSum = c(152,78,163)
RGBcolMolt = c(228,26,28)
alphaLines = .2
alphaPoints = .4
cexPoints = .7
png(filename = paste0(figpath,'trendFits.png'), width = 1020, height = 680) 
  NiMat = M[['Ni']]
  etaVec = M[['eta']]
  par(mar = c(5,5,2,1))
  mnWdiMolt = 0
  for(i in 1:1000) mnWdiMolt = mnWdiMolt + M[['WdiMolt']][[i]][13,]
  mnWdiMolt = mnWdiMolt/1000
  plot(1992:2017 + .6, exp(etaVec[1] + NiMat[,1] + mnWdiMolt), type = 'l', ylim = c(0,12000),
     col = rgb(red = RGBcolRaw[1]/255, green = RGBcolRaw[2]/255,
     blue = RGBcolRaw[3]/255, alpha = alphaLines), 
     cex.axis = 1.5, xlab = 'Year', ylab = 'Counts',
     cex.lab = 2)
  for(i in 2:1000)
    lines(1992:2017 + .6, exp(etaVec[i] + NiMat[,i] + mnWdiMolt), 
			col = rgb(red = RGBcolRaw[1]/255, green = RGBcolRaw[2]/255,
			blue = RGBcolRaw[3]/255, alpha = alphaLines))
  #add sightability
  i = 1
  beta0 = unlist(lapply(MHOaug[['beta']], function(x){x[[1]]}))
  for(i in 1:1000)
    lines(1992:2017 + .6, (exp(etaVec[i] + NiMat[,i] + mnWdiMolt))*
    (1+exp(beta0[i]))/exp(beta0[i]), col = rgb(red = RGBcolEst[1]/255, 
    green = RGBcolEst[2]/255, blue = RGBcolEst[3]/255, alpha = alphaLines))
  x = dgrnd$year + dgrnd$yday/365
  y = dgrnd$adu
  par(mar = c(5,5,2,1))
  points(x[dgrnd$season == 'molt'],y[dgrnd$season == 'molt'],
    pch = 19, cex = cexPoints, col = rgb(red = RGBcolMolt[1]/255, 
    green = RGBcolMolt[2]/255, blue = RGBcolMolt[3]/255, alpha = alphaPoints))
  points(x[dgrnd$season == 'pup'],y[dgrnd$season == 'pup'],
    pch = 19, cex = cexPoints, col = rgb(red = RGBcolAdu[1]/255, 
    green = RGBcolAdu[2]/255, blue = RGBcolAdu[3]/255, alpha = alphaPoints))
  y1 = dgrnd$pup
  points(x[dgrnd$season == 'pup'],y1[dgrnd$season == 'pup'],
    pch = 19, cex = cexPoints, col = rgb(red = RGBcolPup[1]/255, 
    green = RGBcolPup[2]/255, blue = RGBcolPup[3]/255, alpha = alphaPoints))
  points(x[dgrnd$season == 'pup'], y1[dgrnd$season == 'pup'] + 
    y[dgrnd$season == 'pup'], pch = 19, cex = cexPoints,
    col = rgb(red = RGBcolSum[1]/255, 
    green = RGBcolSum[2]/255, blue = RGBcolSum[3]/255, alpha = alphaPoints))

  x2 = dair1$year + dair1$yday/365
  y2 = dair1$aduEst
  points(x2[dair1$season == 'molt'],y2[dair1$season == 'molt'],
    cex = 2, lwd = 2, col = rgb(red = RGBcolMolt[1]/255, 
    green = RGBcolMolt[2]/255, blue = RGBcolMolt[3]/255))
  points(x2[dair1$season == 'pup'],y2[dair1$season == 'pup'],
    col = rgb(red = RGBcolAdu[1]/255, 
    green = RGBcolAdu[2]/255, blue = RGBcolAdu[3]/255), lwd = 2, cex = 2)
  y3 = dair1$pupEst
  points(x2[dair1$season == 'pup'],y3[dair1$season == 'pup'],
    col = rgb(red = RGBcolPup[1]/255, 
    green = RGBcolPup[2]/255, blue = RGBcolPup[3]/255), lwd = 2, cex = 2)
  points(x2[dair1$season == 'pup'],y2[dair1$season == 'pup'] +
    y3[dair1$season == 'pup'], col = rgb(red = RGBcolSum[1]/255, 
    green = RGBcolSum[2]/255, blue = RGBcolSum[3]/255), lwd = 2, cex = 2)
  legend(2006,12000, 
    legend = c('Ground counts, pupping pups',
             'Ground counts, pupping nonpups',
             'Ground counts, pupping sum',
             'Ground counts, molting total',
             'Aerial estimates, pupping pups',
             'Aerial estimates, pupping nonpups',
             'Aerial estimates, pupping sum',
             'Aerial estimates, molting total'),
    pch = c(19,19,19,19,1,1,1,1), 
    col = c(rgb(red = RGBcolPup[1]/255, green = RGBcolPup[2]/255, 
			blue = RGBcolPup[3]/255),
			rgb(red = RGBcolAdu[1]/255, green = RGBcolAdu[2]/255, 
			blue = RGBcolAdu[3]/255),
			rgb(red = RGBcolSum[1]/255, green = RGBcolSum[2]/255, 
			blue = RGBcolSum[3]/255),
			rgb(red = RGBcolMolt[1]/255, green = RGBcolMolt[2]/255, 
			blue = RGBcolMolt[3]/255),
			rgb(red = RGBcolPup[1]/255, green = RGBcolPup[2]/255, 
			blue = RGBcolPup[3]/255),
			rgb(red = RGBcolAdu[1]/255, green = RGBcolAdu[2]/255, 
			blue = RGBcolAdu[3]/255),
			rgb(red = RGBcolSum[1]/255, green = RGBcolSum[2]/255, 
			blue = RGBcolSum[3]/255),
			rgb(red = RGBcolMolt[1]/255, green = RGBcolMolt[2]/255, 
			blue = RGBcolMolt[3]/255)),
    cex = 1.5)

  AbunMat = matrix(0, nrow = 1000, ncol = length(NiMat[,1]))
  for(i in 1:1000)
    AbunMat[i,] = (exp(etaVec[i] + NiMat[,i] + mnWdiMolt))*
      (1+exp(beta0[i]))/exp(beta0[i])
  lines(1992:2017 + .6, apply(AbunMat,2,mean), col = 'black', lwd = 4)
  lines(1992:2017 + .6, apply(AbunMat,2,quantile,probs=.05), 
    col = 'black', lwd = 4, lty = 2)    
  lines(1992:2017 + .6, apply(AbunMat,2,quantile,probs=.95), 
    col = 'black', lwd = 4, lty = 2)  
dev.off()  

## ----fitHOBetaDist, echo=FALSE, include = FALSE, fig.height = 10, fig.width = 10, cache = TRUE----

beta0 = unlist(lapply(MHOaug[['beta']], function(x){x[[1]]}))
png(filename = paste0(figpath,'fitHOBetaDist.png'), width = 680, height = 680) 
  MHO = MHOaug
  par(mar = c(5,5,0,0))
  plot(c(0,1),c(0,45), type = 'n', 
    xlab = 'Proportion', ylab = 'Probability Density',
    cex.axis = 1.5, cex.lab = 2)
  hist(dHOfal$y, add = TRUE, 
    breaks = c(0,.01,.05,.15,.25,.35,.45,.55,.65,.75,.85,.95,.99,1),
    freq = FALSE, col = rgb(1,.1,.0,.2), border = rgb(1,.1,.0,.2))
  makex = c((1:100)/10000,(11:989)/1000,(9900:9999)/10000)
  for(k in 1:length(MHO[[2]])) {
    lines(makex, dbeta(makex,
      exp(beta0[k])/(1+exp(beta0[k]))*MHO[[1]][k],
      (1-exp(beta0[k])/(1+exp(beta0[k])))*MHO[[1]][k]), type = 'l',
      col = rgb(0,0,0,.03))
}
dev.off()

## ----fitHOdate, echo=FALSE, include = FALSE, fig.height = 10, fig.width = 10, cache = TRUE----

beta0 = unlist(lapply(MHOaug[['beta']], function(x){x[[1]]}))
beta1 = unlist(lapply(MHOaug[['beta']], function(x){x[[2]]}))
beta2 = unlist(lapply(MHOaug[['beta']], function(x){x[[3]]}))
png(filename = paste0(figpath,'fitHOdate.png'), width = 680, height = 680) 
  par(mar = c(5,5,0,0))
  plot(c(1,30),c(0,1), type = 'n',
    ylab = 'Haul-Out Probability', xlab = "August Date",
    cex.lab = 2, cex.axis = 1.5)
  x = ((1:30) - 15)/30
  HOvec = matrix(0, nrow = length(beta0), ncol = 30)
  for(k in 1:length(beta0)) {
    Xb = beta0[k] + beta1[k]*x + beta2[k]*x^2
    lines(1:30, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
    HOvec[k,] = exp(Xb)/(1+exp(Xb))
  }
  lines(1:30, apply(HOvec,2,mean), lwd = 4, col = 'blue')
  lines(1:30, apply(HOvec,2,quantile, probs = .05), lty = 2, lwd = 4, col = 'blue')
  lines(1:30, apply(HOvec,2,quantile, probs = .95), lty = 2, lwd = 4, col = 'blue')
dev.off()

## ----postAutoC, echo=FALSE, include = FALSE, fig.height = 10, fig.width = 10, cache = TRUE----

png(filename = paste0(figpath,'postAutoC.png'), width = 750, height = 750)
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  par(mar = c(6,6,6,2))
  plot(density(M[['phiYr']],bw = .01),
    ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
    cex.lab = 2, cex.axis = 1.5, lwd = 3,
    main = 'Yearly, Count Surveys',
    cex.main = 1.7)

  plot(density(M[['phiDyPup']],bw = .02),
    ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
    main = 'Daily, Count Surveys During Pupping', cex.main = 1.7,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)

  plot(density(M[['phiDyMolt']],bw = .05),
    ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
    main = 'Daily, Count Surveys During Molting', cex.main = 1.7,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)

  plot(density(exp(MHOaug[['alpha']])/(1+exp(MHOaug[['alpha']])),bw = .004),
    ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
    main = 'Hourly, Haul-out Model', cex.main = 1.7,
    cex.lab = 2, cex.axis = 1.5, lwd = 3)
dev.off()
