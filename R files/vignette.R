######################################################################
## Compare EpiFilter with APEestim and EpiEstim
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Notes and assumptions
# - load H1N1 influenza 1918 or SARS 2003 data from EpiEstim
# - estimate reproduction numbers with EpiEstim (weekly windows)
# - estimate reproduction numbers with EpiEstim (optimises for prediction)
# - estimate reproduction numbers with EpiFilter (no windows)

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/test/")

# Main packages
library("EpiEstim")
library("caTools")
# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

# Decide if to use pre-processed data (5-day moving average)
preProc = 1; if(preProc){print('Using filtered data')}
# For fixed window get times of incidence curve
win = 7  # weekly

if(!preProc){
  # Load flu data preprocessed with 5-day moving average
  Iflu = read.csv('data/test/Iflu.csv'); Iflu = Iflu[[1]]
  Lflu = read.csv('data/test/Lflu.csv'); Lflu = Lflu[[1]]
  genflu = read.csv('data/test/genflu.csv'); genflu = genflu[[1]]
  
  # Load sars data preprocessed with 5-day moving average
  Isars = read.csv('data/test/Isars.csv'); Isars = Isars[[1]]
  Lsars = read.csv('data/test/Lsars.csv'); Lsars = Lsars[[1]]
  gensars = read.csv('data/test/gensars.csv'); gensars = gensars[[1]]

}else{
  # Load flu data preprocessed with 5-day moving average
  Iflu = read.csv('data/test/IfluFilt.csv'); Iflu = Iflu[[1]]
  Lflu = read.csv('data/test/LfluFilt.csv'); Lflu = Lflu[[1]]
  genflu = read.csv('data/test/genflu.csv'); genflu = genflu[[1]]
  
  # Load sars data preprocessed with 5-day moving average
  Isars = read.csv('data/test/IsarsFilt.csv'); Isars = Isars[[1]]
  Lsars = read.csv('data/test/LsarsFilt.csv'); Lsars = Lsars[[1]]
  gensars = read.csv('data/test/gensars.csv'); gensars = gensars[[1]]
}

# Times of series
tflu = 8:length(Iflu); tsars = 20:length(Isars)

# Plotting of incidence
Ipltflu = Iflu[seq(tflu[1]+1,length(Iflu))]
Ipltsars = Isars[seq(tsars[1]+1,length(Isars))]

######################################################################
## EpiEstim (fixed window) and APEestim (optimised window) computed from:
# Parag KV, Donnelly CA (2020) Using information theory to optimise epidemic models
# for real-time prediction and estimation. PLoS Comput Biol 16(7): e1007990.
######################################################################

# Priors and settings (confidence intervals fixed by a)
Rprior = c(1, 2); a = 0.025
# Clean total infectiousness vectors of NAs
Lflu[is.na(Lflu)] = 0; Lsars[is.na(Lsars)] = 0

# Outputs take list form [[kbest1, modBest1, kbest2, modBest2]]
# Each modBest is also a list of form [[ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci, alpha, beta, pr]]

# Flu results
Rmodflu = apeEstim(Iflu, genflu, Lflu, Rprior, a, tflu[1], 'flu')
# Best estimates and prediction 
plotAPEWindow(Rmodflu[[2]], 'APEestim_flu', Rmodflu[[1]], Ipltflu , folres)

# Specific 7-day window
Rmodflu7 = apeSpecific(Iflu, genflu, Lflu, Rprior, a, tflu[1], win)
plotAPEWindow(Rmodflu7[[2]], 'EpiEstim_flu7', Rmodflu7[[1]], Ipltflu, folres)

# Sars results
Rmodsars = apeEstim(Isars, gensars, Lsars, Rprior, a, tsars[1], 'sars')
# Best estimates and prediction 
plotAPEWindow(Rmodsars[[2]], 'APEestim_sars', Rmodsars[[1]], Ipltsars , folres)

# Specific 7-day window
Rmodsars7 = apeSpecific(Isars, gensars, Lsars, Rprior, a, tsars[1], win)
plotAPEWindow(Rmodsars7[[2]], 'EpiEstim_sars7', Rmodsars7[[1]], Ipltsars, folres)

######################################################################
## EpiFilter: provides formally smoothed and exact estimates
# Method based on Bayesian recursive filtering and smoothing
######################################################################

# Setup grid and noise parameters
Rmin = 0.01; Rmax = 10; eta = 0.1

# Uniform prior over grid of size m
m = 200; pR0 = (1/m)*rep(1, m)
# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)

# Time series lengths
nflu = length(tflu); nsars = length(tsars)

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt_flu = epiFilter(Rgrid, m, eta, pR0, nflu, Lflu[tflu], Iflu[tflu], a)
Rfilt_sars = epiFilter(Rgrid, m, eta, pR0, nsars, Lsars[tsars], Isars[tsars], a)

# Causal predictions from filtered estimates [pred predci]
Ifilt_flu = recursPredict(Rgrid, Rfilt_flu[[4]], Lflu[tflu], Rfilt_flu[[3]], a)
Ifilt_sars = recursPredict(Rgrid, Rfilt_sars[[4]], Lsars[tsars], Rfilt_sars[[3]], a)

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth_flu = epiSmoother(Rgrid, m, Rfilt_flu[[4]], Rfilt_flu[[5]], nflu, Rfilt_flu[[6]], a)
Rsmooth_sars = epiSmoother(Rgrid, m, Rfilt_sars[[4]], Rfilt_sars[[5]], nsars, Rfilt_sars[[6]], a)

# Smoothed predictions from filtered estimates [pred predci]
Ismooth_flu = recursPredict(Rgrid, Rsmooth_flu[[4]], Lflu[tflu], Rsmooth_flu[[3]], a)
Ismooth_sars = recursPredict(Rgrid, Rsmooth_sars[[4]], Lsars[tsars], Rsmooth_sars[[3]], a)

 # Plot estimates and predictions from filtering
plotEpiFilter(Rfilt_flu[[3]][2:nflu], Rfilt_flu[[2]][, 2:nflu], Ifilt_flu[[1]], Ifilt_flu[[2]],
              'EpiFilter_flu', Ipltflu, folres, eta)
plotEpiFilter(Rfilt_sars[[3]][2:nsars], Rfilt_sars[[2]][, 2:nsars], Ifilt_sars[[1]], Ifilt_sars[[2]],
              'EpiFilter_sars', Ipltsars, folres, eta)

# Plot estimates and predictions from smoothing
plotEpiFilter(Rsmooth_flu[[3]][2:nflu], Rsmooth_flu[[2]][, 2:nflu], Ismooth_flu[[1]], Ismooth_flu[[2]],
              'EpiSmooth_flu', Ipltflu, folres, eta)
plotEpiFilter(Rsmooth_sars[[3]][2:nsars], Rsmooth_sars[[2]][, 2:nsars], Ismooth_sars[[1]], Ismooth_sars[[2]],
              'EpiSmooth_sars', Ipltsars, folres, eta)

