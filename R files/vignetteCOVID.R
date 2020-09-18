######################################################################
## Apply EpiFilter to COVID data from New Zealand
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Notes and assumptions
# - load COVID incidence curves from WHO dashboard
# - estimate reproduction numbers with EpiFilter 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/covid/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

# Load data from WHO for New Zealand
alldata = read.csv('data/covid/WHO-COVID-19-global-data.csv')
idcountry = which(alldata$Country == 'New Zealand')

# Incidence and dates
Iday = alldata$New_cases[idcountry]
dates  = alldata$Date_reported[idcountry]
# Time series lengths
nday = length(dates); tday = 1:nday

# Approxumate serial interval distribution from Ferguson et al
wdist = dgamma(tday, shape = 2.3669, scale = 2.7463)

# Total infectiousness
Lday = rep(0, nday) 
for(i in 2:nday){
  # Total infectiousness
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}

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

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

 # Plot estimates and predictions from filtering
plotEpiFilter(Rfilt[[3]][2:nday], Rfilt[[2]][, 2:nday], Ifilt[[1]], Ifilt[[2]],
              'EpiFilter', Iday[2:nday], folres, eta)

# Plot estimates and predictions from smoothing
plotEpiFilter(Rsmooth[[3]][2:nday], Rsmooth[[2]][, 2:nday], Ismooth[[1]], Ismooth[[2]],
              'EpiSmooth', Iday[2:nday], folres, eta)

