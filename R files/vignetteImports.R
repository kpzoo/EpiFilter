######################################################################
## Local (R, Z) framework for investigating transmission in New Zealand
# R is smoothed local reprod. number i.e. we account for imports/local cases
# Z is confidence in elimination in %, also accounting for local/imports

# Framework merges methods from:
# Parag K. (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves”, BioRxiv 2020.09.14.20194589.
# Parag K, Donnelly C, et al. (2020) “An exact method for quantifying the reliability
# of end-of-epidemic declarations in real time”, PLOS Comput Biol, In Press.

######################################################################

# Notes and assumptions
# - COVID incidence curves (local/imported) from EpiSurv
# - estimate reproduction numbers with EpiFilter 
# - compute associated confidence in end-of-epidemic declarations
# - do not start incidence curve with 0 total cases

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# For dates from excel
library(lubridate); library(openxlsx)
# For plotting area charts
library(plotly)

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/imports/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

######################################################################
## New Zealand data cleaning and renewal model setup
######################################################################

# Load data from EpiSurv for New Zealand
alldata = read.xlsx('data/imports/covid-cases-07oct20_confirmed.xlsx', detectDates = TRUE)
lendata = length(alldata$`Confirmed.Covid-19.cases`)

# Case dates (of reporting) after truncating header
dates = alldata$`Confirmed.Covid-19.cases`[3:lendata]
# Range of all possible dates
datemin = min(dates); datemax = max(dates)
# All dates that form incidence curves
tdates = seq(as.Date(datemin), as.Date(datemax), by = "day")
nday = length(tdates); allcases = lendata-2

# Corresponding travel history of these cases
travelled = alldata$X5[3:lendata]
# Classify into imported vs local cases
local = dates[travelled == "No"]; len1 = length(local)
imported = dates[travelled == "Yes"]; len2 = length(imported)
# Cases with unknown travel history
unknown = dates[travelled != "No" & travelled != "Yes"]

# Construct local and imported incidence curves
Iloc = rep(0, nday); Iintro = Iloc
for (i in 1:nday) {
  Iloc[i] = length(which(local == tdates[i]))
  Iintro[i] = length(which(imported == tdates[i]))
}

# Total cases for computing total infectiousness
Iday = Iloc + Iintro
# Check accounted for all cases
if(allcases != sum(Iday) + length(unknown)){
  stop("All cases not accounted for")
}

# Time series indices and shape-scale parameters
tday = 1:nday; pms = c(2.3669, 2.7463)
# Approximate serial interval distribution from Ferguson et al
tdist = c(0, tday); wdist = rep(0, nday)
for (i in 1:nday){
  wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) - 
    pgamma(tdist[i], shape = pms[1], scale = pms[2])
}

# Total infectiousness
Lday = rep(0, nday) 
for(i in 2:nday){
  # Total infectiousness
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}

# Local and imported cases
pdf(file=paste0(folres, 'incidence', '.pdf')) 
plot(tdates, Iday, type = 'h', lwd = 2, col='grey', bty = 'l',
     xlab = 'time (days)', ylab = 'incidence')
lines(tdates, Iloc, type = 'h', col = 'red')
dev.off()

######################################################################
## EpiFilter: provides formally smoothed and exact estimates
# Method based on Bayesian recursive filtering and smoothing
######################################################################

# Setup grid and noise parameters
Rmin = 0.01; Rmax = 10; eta = 0.1

# Uniform prior over grid of size m
m = 500; pR0 = (1/m)*rep(1, m)
# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iloc[tday], 0.025)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredictAdj(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025, 2*max(Iloc))

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredictAdj(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025, 2*max(Iloc))

 # Plot estimates and predictions from filtering
plotEpiFilter(Rfilt[[3]][2:nday], Rfilt[[2]][, 2:nday], Ifilt[[1]], Ifilt[[2]],
              'EpiFilter', Iloc[2:nday], folres, eta)
# Plot estimates and predictions from smoothing
plotEpiFilter(Rsmooth[[3]][2:nday], Rsmooth[[2]][, 2:nday], Ismooth[[1]], Ismooth[[2]],
              'EpiSmooth', Iloc[2:nday], folres, eta)

######################################################################
## Compute resulting Z number i.e confidence in the end of the epidemic
######################################################################

# Filtered and smoothed Z numbers
tZ = 2:length(Iday); lenZ = length(tZ)
ZnumS = rep(0, lenZ); ZnumF = ZnumS

# Compute Z scores iteratively
for(i in 1:lenZ){
  # Main elimination probability function
  Zout = computeZ(Iloc[1:tZ[i]], Iday[1:tZ[i]], Rgrid, m, pms, eta, pR0)
  ZnumF[i] = Zout[1]; ZnumS[i] = Zout[2]
}

# Z scores with time
pdf(file=paste0(folres, 'znumbers', '.pdf')) 
plot(tdates[2:length(Iday)], ZnumF, type = 'l', lwd = 2, col='grey', bty = 'l',
     xlab = 'time (days)', ylab = '% confidence in elimination')
lines(tdates[2:length(Iday)], ZnumS, type = 'l', col = 'red', lwd = 2)
dev.off()

# Dates of 95% and 99% declaration
id95 = c(which(ZnumF >= 95)[1], which(ZnumS >= 95)[1])
id99 = c(which(ZnumF >= 99)[1], which(ZnumS >= 99)[1])
# Add 1 as ids referenced to 2:length(Iday)
id95 = id95 + 1; id99 = id99 + 1
t95 = tdates[id95]; t99 = tdates[id99]

######################################################################
## Alignment of policy with transmission and elimination
######################################################################

# Time of lockdown enforcement (level 4 alert)
tlock = "2020-03-25"; idlock = which(tdates == tlock)
# Time when NPIs are relaxed (alert level 2)
trelax = "2020-05-14"; idrelax = which(tdates == trelax)
# Time that first wave id declared eliminated (alert level 1)
telim = "2020-06-09"; idelim = which(tdates == telim)
# Time when NPIs enforced against resurgence (alert level 3)
tsecwv = "2020-08-13"; idsecwv = which(tdates == tsecwv)
# Last analysed time (NZ believes COVID is under control again)
tpres = "2020-10-07"; idpres = which(tdates == tpres)

# Plot of timeline of incidence
ymax = max(max(Iloc), max(Iintro)) + 5
# Local and imported cases (not stacked)
fig = plot_ly(x = tdates, y = Iloc, type = 'scatter', mode = 'lines', name = 'local', fill = 'tozeroy')
fig = fig %>% add_trace(x = tdates, y = Iintro, name = 'imported', fill = 'tozeroy')
# All policy dates
fig = fig %>% add_trace(x = c(tlock, tlock), y = c(0, ymax), name = 'lockdown',
                        showlegend = FALSE, line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(trelax, trelax), y = c(0, ymax), name = 'relaxed NPIs',
                        showlegend = FALSE, line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(telim, telim), y = c(0, ymax), name = 'elimination',
                        showlegend = FALSE, line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(tsecwv, tsecwv), y = c(0, ymax), name = 'enforced NPIs',
                        showlegend = FALSE, line = list(color = 'grey'))
fig = fig %>% layout(xaxis = list(title = 'time (days)'), yaxis = list(title = 'incidence'))
figIncid = fig;

# R estimates with confidence limits
ymax = max(Rsmooth[[2]][2,]) + 0.2
# Smoothed local R estimates (95% confidence intervals)
fig = plot_ly(x = tdates, y = Rsmooth[[2]][1,], type = 'scatter', mode = 'lines', name = 'lower CI (R)',
              showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = Rsmooth[[2]][2,], name = 'upper CI (R)', fill = 'tonexty',
                        showlegend = FALSE, line = list(color = 'red'))
fig = fig %>% add_trace(x = tdates, y = rep(1, nday), name = 'R = 1', line = list(color = 'lightblue'))
# All policy dates
fig = fig %>% add_trace(x = c(tlock, tlock), y = c(0, ymax), name = 'lockdown', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(trelax, trelax), y = c(0, ymax), name = 'relaxed NPIs', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(telim, telim), y = c(0, ymax), name = 'elimination', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(tsecwv, tsecwv), y = c(0, ymax), name = 'enforced NPIs', line = list(color = 'grey'))
fig = fig %>% layout(xaxis = list(title = 'time (days)' ), yaxis = list(title = 'local reprod. number'))
figRest = fig;

# Z estimates 
fig = plot_ly(x = tdates, y = c(0, ZnumS), type = 'scatter', mode = 'lines', name = 'Z number',
              showlegend = FALSE, line = list(color = 'dodgerblue'))
fig = fig %>% add_trace(x = tdates, y = rep(95, nday), name = 'Z = 95', line = list(color = 'lightblue'))
# All policy dates
ymax = 105
fig = fig %>% add_trace(x = c(tlock, tlock), y = c(0, ymax), name = 'lockdown', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(trelax, trelax), y = c(0, ymax), name = 'relaxed NPIs', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(telim, telim), y = c(0, ymax), name = 'elimination', line = list(color = 'grey'))
fig = fig %>% add_trace(x = c(tsecwv, tsecwv), y = c(0, ymax), name = 'enforced NPIs', line = list(color = 'grey'))
fig = fig %>% layout(xaxis = list(title = 'time (days)' ), yaxis = list(title = '% confidence in elimination'))
figZest = fig;


