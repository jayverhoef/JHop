library(akpvdata)
library(akpvana)


########################################################################
# -------------------- HAULOUT GRAPHICS --------------------------------
########################################################################

  data(sHOe_GLAC)
  MHO = sHOe_GLAC_spr$M
  dHO = sHOe_GLAC_spr$dHO
  
 
# ------------- MCMC traces for regression parameters ------------------

  X11()
  layout(matrix(1:7, ncol = 1), heights = c(1.09,1,1,1,1,1,1.48))
  par(mar = c(0,5,1,1))
  beta0 = unlist(lapply(MHO[['beta']], function(x){x[[1]]}))
  plot(beta0, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[0]), xaxt = 'n')
  beta1 = unlist(lapply(MHO[['beta']], function(x){x[[2]]}))
  par(mar = c(0,5,0,1))
  plot(beta1, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[1]), xaxt = 'n')
  beta2 = unlist(lapply(MHO[['beta']], function(x){x[[3]]}))
  plot(beta2, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[2]), xaxt = 'n')
  beta3 = unlist(lapply(MHO[['beta']], function(x){x[[4]]}))
  plot(beta3, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[3]), xaxt = 'n')
  beta4 = unlist(lapply(MHO[['beta']], function(x){x[[5]]}))
  plot(beta4, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[4]), xaxt = 'n')
  beta5 = unlist(lapply(MHO[['beta']], function(x){x[[6]]}))
  plot(beta5, type = 'l', xlab = '', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[5]), xaxt = 'n')
  beta6 = unlist(lapply(MHO[['beta']], function(x){x[[7]]}))
  par(mar = c(5,5,0,1))
  plot(beta6, type = 'l', xlab = 'MCMC Index', cex.lab = 2, cex.axis = 1.5,
    ylab = expression(beta[6]))

# ----------Last realization of the temporally-autocorrelated  ---------
# --------------------random effects per seal---------------------------
# -----------Compared to observed haul-out proportion-------------------

  #names and sample sizes
  table(dHO$speno)

debug(HOautocorrRErealizeFit)
  HOautocorrRErealizeFit(MHO[['eps']], dHO, 1)
  HOautocorrRErealizeFit(MHO[['eps']], dHO, 2) #etc.


# --------------------------- Date effect ------------------------------

X11()
par(mar = c(5,5,0,0))
plot(c(1,75),c(0,.8), type = 'n',
  ylab = 'Haul-Out Probability', xlab = "Days since 15 July",
  cex.lab = 2, cex.axis = 1.5)
x = ((1:75) - 30)/30
HOvec = matrix(0, nrow = length(beta1), ncol = 75)
for(k in 1:length(beta1)) {
  Xb = beta0[k] + beta1[k]*x + beta2[k]*x^2
  lines(1:75, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
  HOvec[k,] = exp(Xb)/(1+exp(Xb))
}
lines(1:75, apply(HOvec,2,mean), lwd = 4, col = 'blue')
lines(1:75, apply(HOvec,2,quantile, probs = .05), lty = 2, lwd = 4, 
  col = 'blue')
lines(1:75, apply(HOvec,2,quantile, probs = .95), lty = 2, lwd = 4, 
  col = 'blue')
hist((dHO$dy-as.POSIXlt("2004-05-01")$yday), add = TRUE, freq = FALSE,
	col = rgb(.8,.1,.1,.4))
hist(dstk$day-as.POSIXlt(as.POSIXct("2005-07-15"))$yday, freq = FALSE, 
	add = TRUE, col = rgb(.1,.8,.1,.4))

# graphics - tide effect
X11()
par(mar = c(5,5,5,0))
plot(c(-5,5),c(0,.8), type = 'n',
  ylab = 'Haul-Out Probability', xlab = "Hours from Low Tide",
  cex.lab = 2, cex.axis = 1.5)
x = (-50:50)/25
HOvec = matrix(0, nrow = length(beta3), ncol = length(x))
for(k in 1:length(beta3)) {
  Xb = beta0[k] + beta3[k]*x + beta4[k]*x^2
  lines((-50:50)/10, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
  HOvec[k,] = exp(Xb)/(1+exp(Xb))
}
lines((-50:50)/10, apply(HOvec,2,mean), lwd = 4, col = 'blue')
lines((-50:50)/10, apply(HOvec,2,quantile, probs = .05), lty = 2, lwd = 4, col = 'blue')
lines((-50:50)/10, apply(HOvec,2,quantile, probs = .95), lty = 2, lwd = 4, col = 'blue')
hist(dHO$minutes_from_low/60, add = TRUE, freq = FALSE,
	col = rgb(.8,.1,.1,.4))
hist(dstk$time_from_low/60, add = TRUE, freq = FALSE,
	col = rgb(.1,.8,.1,.4))

# graphics - hour-of-day effect
X11()
par(mar = c(5,5,0,0))
plot(c(-12,12),c(0,.8), type = 'n',
  ylab = 'Haul-Out Probability', xlab = "Hours from Solar Noon",
  cex.lab = 2, cex.axis = 1.5)
x = (-60:60)/30
HOvec = matrix(0, nrow = length(beta1), ncol = length(x))
for(k in 1:length(beta1)) {
  Xb = beta0[k] + beta5[k]*x + beta6[k]*x^2
  lines((-60:60)/5, exp(Xb)/(1+exp(Xb)), col = rgb(0,0,0,.03))
  HOvec[k,] = exp(Xb)/(1+exp(Xb))
}
lines((-60:60)/5, apply(HOvec,2,mean), lwd = 4, col = 'blue')
lines((-60:60)/5, apply(HOvec,2,quantile, probs = .05), lty = 2, lwd = 4, col = 'blue')
lines((-60:60)/5, apply(HOvec,2,quantile, probs = .95), lty = 2, lwd = 4, col = 'blue')
hist(dHO$solhr - 12, add = TRUE, freq = FALSE,
	col = rgb(.8,.1,.1,.4))
hist(dstk$hr - 12, add = TRUE, freq = FALSE,
	col = rgb(.1,.8,.1,.4))

plot(c(0,1),c(0,45), type = 'n', 
	xlab = 'Proportion', ylab = 'Probability Density',
	cex.lab = 2, cex.axis = 1.5)
hist(dHO$y, add = TRUE, 
	breaks = c(0,.01,.05,.15,.25,.35,.45,.55,.65,.75,.85,.95,.99,1),
	freq = FALSE, col = rgb(1,.1,.0,.2), border = rgb(1,.1,.0,.2))
makex = c((1:100)/10000,(11:989)/1000,(9900:9999)/10000)
for(k in 1:length(M[[2]])) {
lines(makex, dbeta(makex,
	exp(beta0[k])/(1+exp(beta0[k]))*M[[1]][k],
	(1-exp(beta0[k])/(1+exp(beta0[k])))*M[[1]][k]), type = 'l',
	col = rgb(0,0,0,.03))
}

# trace of ith animal's random effect
i = 4
gam0 = unlist(lapply(MHO[['gam']], function(x){x[[i]]}))
plot(gam0, type = 'l')

#all animal's random effects as posterior densities
par(mar = c(5,5,.1,.1))
plot(c(-2.5,2.5),c(0,4), type = 'n', main = '',
	xlab = 'Random Effect Value', ylab = 'Posterior Density',
	cex.lab = 2, cex.axis = 1.5)
for(i in 1:length(MHO[['gam']][[1]])) {
	gam0 = unlist(lapply(MHO[['gam']], function(x){x[[i]]}))
	lines(density(gam0, bw = .08), col = rgb(0,0,1,.4), lwd = 3)
}

#posterior density of hourly Haul-out autocorrelation
X11()
par(mar = c(5,5,5,1))
plot(density(exp(MHO[['alpha']])/(1+exp(MHO[['alpha']])),bw = .01),
  ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(MHO[['alpha']], type = 'l')

#posterior density of variance parameters
par(mar = c(5,5,5,1))
plot(density(M[['sigGam']],bw = .1),
  ylab = 'Posterior Density', xlab = 'Seal Random Effect Variance',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(MHO[['sigGam']], type = 'l')

par(mar = c(5,5,5,1))
plot(density(M[['sigmaARexp']],bw = .01),
  ylab = 'Posterior Density', xlab = 'Haul-out Autocorrelation Variance',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(M[['sigmaARexp']], type = 'l')

X11()
par(mar = c(5,5,5,1))
plot(density(M[['sigmaAR1']],bw = .01),
  ylab = 'Posterior Density', xlab = 'Autocorrelated Error Variance',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(M[['sigmaAR1']], type = 'l')

X11()
par(mar = c(5,5,5,1))
plot(density(M[['sigmaInd']],bw = .01),
  ylab = 'Posterior Density', xlab = 'Independent Error Variance',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(M[['sigmaInde']], type = 'l')

par(mar = c(5,5,5,1))
plot(density(M[['sigEta']],bw = .01),
  ylab = 'Posterior Density', xlab = 'Count Independence Variance',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(M[['sigEta']], type = 'l')


#posterior density of annual count autocorrelation
X11()
par(mar = c(5,5,5,1))
plot(density(exp(-M[['rho']]),bw = .02),
  ylab = 'Posterior Density', xlab = 'Lag-1 Autocorrelation',
  main = 'Hourly, Haul-out Model', cex.main = 2,
  cex.lab = 2, cex.axis = 1.5, lwd = 3)
plot(M[['rho']], type = 'l')


# 1000 corrections factors for 15 August, low tide, solar noon
X11()
cf = (1+exp(beta0))/exp(beta0)
plot(density(cf,bw = .05), main = '',
  ylab = 'Posterior Density', xlab = 'Standardized Correction Factor',
  cex.lab = 2, cex.axis = 1.5, lwd = 3)

#matrix of all tau values
X11()
tauMat = NULL
for(i in 1:length(M[['tau']]))  
	tauMat = cbind(tauMat,M[['tau']][[i]])
# trace of ith sites random effect
sqrt(apply(tauMat,1,var))
i = 8
plot(tauMat[i,], type = 'l')
plot(exp(tauMat[i,]), type = 'l')
plot(exp(tauMat[i,1000] + M[['del']][[1000]][i,]), pch = 19, type = 'l', col = 'blue')

# N condition on site, by year and MCMC
X11()
NiMat = NULL
# jth site
j = 5
for(i in 1:length(M[['Nmat']]))  
	NiMat = rbind(NiMat,M[['Nmat']][[i]][j,])
# trace of kth year
k = 7
plot(NiMat[,k], type = 'l')
apply(NiMat,1,mean)

dstk[dstk$polyid == levels(dstk$polyid)[8],'count']
aggregate(dstk$count,list(dstk$polyid),mean)

# N condition on year, by site and MCMC
X11()
NjMat = NULL
# jth year
j = 9
for(i in 1:length(M[['Nmat']]))  
	NjMat = rbind(NjMat,M[['Nmat']][[i]][,j])
# trace of kth site
k = 14
plot(NjMat[,k], type = 'l')
apply(NjMat,2,mean)

# realization of 3 delta values for ith iteration
X11()
site = 14
del1 = M[['del']][[1]]
del2 = M[['del']][[50]]
del3 = M[['del']][[100]]
plot(c(1,20), c(min(del1,del2,del3), max(del1,del2,del3)), type = 'n')
lines(1:20, del1[site,], col = rgb(.8,.1,.1))
points(1:20, del1[site,], pch = 19, col = rgb(.8,.1,.1))
lines(1:20, del2[site,], col = rgb(.1,.8,.1))
points(1:20, del2[site,], pch = 19, col = rgb(.1,.8,.1))
lines(1:20, del3[site,], col = rgb(.1,.1,.8))
points(1:20, del3[site,], pch = 19, col = rgb(.1,.1,.8))

var(as.vector(M[['Emat']][[100]]))
plot(M[['sigEta']], type = 'l')

X11()
pop = NULL
for(i in 1:length(M[['Nmat']])){
junk = M[['Nmat']][[i]]
	pop = rbind(pop,apply(junk,2,sum))
}
bot = apply(pop,2,quantile, prob = .025)
top = apply(pop,2,quantile, prob = .975)
plot(c(1,dim(M[['del']][[1]])[2]), c(min(bot),max(top)), type = 'n',
	xaxt = 'n', ylab = 'Estimated Standardized Count',
	xlab = '')
axis(1, at = c(5,10,15,20), labels = c(2000, 2005, 2010, 2015))
points(apply(pop,2,mean), pch = 19, cex = 2)
for(i in 1:dim(M[['del']][[1]])[2])
	lines(c(i,i),c(bot[i],top[i]), lty = 1, lwd = 2)

table(dstk$yr)


aggregate(dstk$count,list(dstk$polyid),function(v){round(mean(v),1)})


png('/media/jay/ExtraDrive1/00NMML/activePapers/HSsurv2018/HSsurv2018_package/HSsurv2018/inst/scripts/sitebysite1-12.png',
  width = 960, height = 1440)
layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))
for(site in 1:12) {
# condition on site, matrix of year and MCMC value
pop = NULL
cnts = dstk[dstk$polyid == levels(dstk$polyid)[site],c('yr','count')]
for(i in 1:length(M[['Nmat']])){
junk = M[['Nmat']][[i]][site,]
	pop = rbind(pop, junk)
}
par(mar = c(5,5,5,1))
plot(c(1,dim(M[['Nmat']][[1]])[2]), c(min(pop,cnts[,2]),max(pop,cnts[,2])), type = 'n',
	xaxt = 'n', ylab = 'Abundance Estimate', cex.main = 2,
	xlab = '', main = paste('site', site), cex.lab = 2, cex.axis = 1.5)
axis(1, at = c(5,10,15,20), labels = c(2000, 2005, 2010, 2015))
for(i in 1:length(M[['Nmat']]))
	lines(1:dim(M[['Nmat']][[1]])[2],pop[i,], lty = 1, lwd = 2, col = rgb(0,0,0,.03))
points(cnts[,1]-1995, cnts[,2], pch = 19, cex = 3, col = 'blue')
}
dev.off()

png('/media/jay/ExtraDrive1/00NMML/activePapers/HSsurv2018/HSsurv2018_package/HSsurv2018/inst/scripts/sitebysite13-24.png',
  width = 960, height = 1440)
layout(matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE))
for(site in 13:24) {
# condition on site, matrix of year and MCMC value
pop = NULL
cnts = dstk[dstk$polyid == levels(dstk$polyid)[site],c('yr','count')]
for(i in 1:length(M[['Nmat']])){
junk = M[['Nmat']][[i]][site,]
	pop = rbind(pop, junk)
}
par(mar = c(5,5,5,1))
plot(c(1,dim(M[['Nmat']][[1]])[2]), c(min(pop,cnts[,2]),max(pop,cnts[,2])), type = 'n',
	xaxt = 'n', ylab = 'Abundance Estimate', cex.lab = 2, cex.axis = 1.5,
	xlab = '', main = paste('site', site), cex.main = 2)
axis(1, at = c(5,10,15,20), labels = c(2000, 2005, 2010, 2015))
for(i in 1:length(M[['Nmat']]))
	lines(1:dim(M[['Nmat']][[1]])[2],pop[i,], lty = 1, lwd = 2, col = rgb(0,0,0,.03))
points(cnts[,1]-1995, cnts[,2], pch = 19, cex = 3, col = 'blue')
}
dev.off()

png('/media/jay/ExtraDrive1/00NMML/activePapers/HSsurv2018/HSsurv2018_package/HSsurv2018/inst/scripts/sitebysite4DenMod.png',
  width = 960, height = 640)
layout(matrix(1:12, nrow = 2, ncol = 2, byrow = TRUE))
for(site in c(1,8,14,23)) {
# condition on site, matrix of year and MCMC value
pop = NULL
cnts = dstk[dstk$polyid == levels(dstk$polyid)[site],c('yr','count')]
for(i in 1:length(M[['Nmat']])){
junk = M[['Nmat']][[i]][site,]
	pop = rbind(pop, junk)
}
par(mar = c(5,5,5,1))
plot(c(1,dim(M[['Nmat']][[1]])[2]), c(min(pop,cnts[,2]),max(pop,cnts[,2])), type = 'n',
	xaxt = 'n', ylab = 'Abundance Estimate', cex.lab = 2, cex.axis = 1.5,
	xlab = '', main = paste('site', site), cex.main = 2)
axis(1, at = c(5,10,15,20), labels = c(2000, 2005, 2010, 2015))
for(i in 1:length(M[['Nmat']]))
	lines(1:dim(M[['Nmat']][[1]])[2],pop[i,], lty = 1, lwd = 2, col = rgb(0,0,0,.03))
points(cnts[,1]-1995, cnts[,2], pch = 19, cex = 3, col = 'blue')
}
dev.off()

dstk[dstk$polyid == levels(dstk$polyid)[site],c('yr','count','daystd','hrstd','tide')]

pop1 = 0*del
for(i in 1:length(M[['Nmat']])){
junk = M[['Nmat']][[i]]
	pop1 = pop1 + junk
}
pop1 = pop1/length(M[['Nmat']])
apply(pop1,2,sum)

aggregate(dstk$count,list(dstk$polyid),mean)

X11()
D = outer(1:dim(M[['del']][[1]])[2],rep(1,times = dim(M[['del']][[1]])[2]))
D = abs(D - t(D))
site = 8
del1 = M[['del']][[1]]
del2 = M[['del']][[500]]
del3 = M[['del']][[1000]]
L = eigen(exp(-(D*M[['rho']][1])) + diag(rep(1e-10, times = dim(M[['del']][[1]])[2])))
L1 = L$vectors
L = eigen(exp(-(D*M[['rho']][500])) + diag(rep(1e-10, times = dim(M[['del']][[1]])[2])))
L2 = L$vectors
L = eigen(exp(-(D*M[['rho']][1000])) + diag(rep(1e-10, times = dim(M[['del']][[1]])[2])))
L3 = L$vectors
plot(c(1,dim(M[['del']][[1]])[2]), c(min(L1 %*% del1[site,],L2 %*% del2[site,], 
  L3 %*% del3[site,]), 
  max(L1 %*% del1[site,],L2 %*% del2[site,], 
  L3 %*% del3[site,])), type = 'n', ylab = 'Smooth Temporal Fit',
  xlab = 'Years from 1995')
lines(1:dim(M[['del']][[1]])[2], L1 %*% del1[site,], col = rgb(.8,.1,.1))
points(1:dim(M[['del']][[1]])[2], L1 %*% del1[site,], pch = 19, col = rgb(.8,.1,.1))
lines(1:dim(M[['del']][[1]])[2], L2 %*% del2[site,], col = rgb(.1,.8,.1))
points(1:dim(M[['del']][[1]])[2], L2 %*% del2[site,], pch = 19, col = rgb(.1,.8,.1))
lines(1:dim(M[['del']][[1]])[2], L3 %*% del3[site,], col = rgb(.1,.1,.8))
points(1:dim(M[['del']][[1]])[2], L3 %*% del3[site,], pch = 19, col = rgb(.1,.1,.8))


site = 8
rho = M[['rho']][length(M[['rho']])]
sigmaAR1 = M[['sigmaAR1']][length(M[['rho']])]
sigmaInde = M[['sigmaInde']][[length(M[['rho']])]]
sigmaREday = M[['sigmaREday']][[length(M[['rho']])]]
Nmat = M[['Nmat']][[length(M[['Nmat']])]]
REday = M[['REday']][[length(M[['REday']])]]
X11()
plot(Nmat[site,], type = 'l', lwd = 3)
del = M[['del']][[length(M[['del']])]]
inde = M[['inde']][[length(M[['inde']])]]
lines(exp(inde[site,]), col = 'blue', lwd = 3)
tau = M[['tau']][[length(M[['tau']])]]
lines(exp(tau[site] + del[site,]), col = 'green2', lwd = 3)
X11()
plot(M[['rho']], type = 'l')
X11()
plot(M[['sigmaAR1']], type = 'l')
X11()
plot(M[['sigmaREday']], type = 'l')
X11()
plot(unlist(lapply(M[['sigmaInde']], function(x){x[[site]]})), type = 'l')
X11()
plot(unlist(lapply(M[['tau']], function(x){x[[site]]})), type = 'l')
X11()
plot(exp(unlist(lapply(M[['tau']], function(x){x[[site]]}))), type = 'l')
X11()
plot(unlist(lapply(M[['REday']], function(x){x[[2]]})), type = 'l')


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

X11()
plot(dpois(range, exp(tau[i] + (L1 %*% del[i,1:nevecs])[j]),log = TRUE)) 
X11()
plot(dbinom(cnt, range, p_i, log = TRUE))
X11()
plot(dpois(range, exp(tau[i] + (L1 %*% del[i,1:nevecs])[j]),log = TRUE) +
dbinom(cnt, range, p_i, log = TRUE))

X11()
plot(dpois(range, exp(tau[i] + (L1 %*% del[i,1:nevecs])[j]),log = TRUE)) 
X11()
plot(dbebivhat(as.integer(cnt), range, p_i, vht, log = TRUE))
X11()
plot(dpois(range, exp(tau[i] + (L1 %*% del[i,1:nevecs])[j]),log = TRUE) + 
 dbebivhat(as.integer(cnt), range, p_i, vht, log = TRUE))
