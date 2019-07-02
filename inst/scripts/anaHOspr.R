library(BBPa)
library(JHop)
basedir1 = '/media/jay/Hitachi2GB/00NMML/ActiveRPack/'
basedir2 = 'Harb_Seal_Surv/JHop/data/'
basedir = paste0(basedir1, basedir2)


#-----------------------------------------------------------------------
# Glacial
#-----------------------------------------------------------------------
#initialize parameters with list of stocks being grouped
hi = haulout_initial_JHop('spr')
#fit model with MCMC
sHOe_GLAC_spr = stock_haulout_estimation(hi, she = NULL, ninit = 1000, 
  nburnloops = 100, nperburn = 100, niter = 10000, thin = 100)
#save to analysis package for further analysis
save(sHOe_GLAC_spr, file = paste0(basedir,'sHOe_GLAC_spr.rda'))

#can start MCMC from previous output if desired
#data(sHOe_GLAC)
#hi = haulout_initial('glac')
sHOe_GLAC_spr = stock_haulout_estimation(hi, she = sHOe_GLAC_spr,
  ninit = 10, nburnloops = 20, nperburn = 100, niter = 100000, thin = 100)
#save to analysis package for further analysis
save(sHOe_GLAC_spr, file = paste0(basedir,'sHOe_GLAC_spr.rda'))
