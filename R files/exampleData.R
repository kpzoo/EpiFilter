######################################################################
## Demonstrate APE with EpiEstim package
######################################################################

# Assumptions and modifications
# - compute EpiEstim for influenza and SARS data
# - get APE-regularised estimates for comparison

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Main packages
library("EpiEstim")
library("caTools")

# Key functions
source('apeEstim.R')
source('apePredPost.R')
source('apeSpecific.R')
source('plotAPEWindow.R')

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/")

# Load data on pandemic flu in Baltimore in 1918
data("Flu1918");
Iflu = Flu1918$incidence
genflu = Flu1918$si_distr
# Total infectiousness
Lflu = overall_infectivity(Iflu, genflu)

# Load data on SARS in Hong Kong in 2003
data("SARS2003")
Isars = SARS2003$incidence
gensars = SARS2003$si_distr
# Total infectiousness
Lsars = overall_infectivity(Isars, gensars)

######################################################################
## Conventional EpiEstim estimates
######################################################################

# Start and end times for weekly window in flu
t_start = seq(2, length(Iflu)-6); t_end = t_start + 6
# Estimates on a 7 day window - considered best in Cori 2013 for this data
estflu = estimate_R(Iflu, method ="non_parametric_si", config = make_config(list(
  si_distr = genflu, t_start = t_start, t_end = t_end)))

# Extract outputs
tflu = estflu$R$t_end # end of window
Rflu = estflu$R$`Mean(R)`
RfluCI = matrix(NA, nrow = 2, ncol = length(Rflu))
RfluCI[1,] = estflu$R$`Quantile.0.025(R)`; RfluCI[2,] = estflu$R$`Quantile.0.975(R)`

# Start and end times for weekly window in SARS
t_start = seq(14, length(Isars)-6)
t_end = t_start + 6
# Estimates on a 7 day window - considered best in Cori 2013 for this data
estsars = estimate_R(Isars, method ="non_parametric_si", config = make_config(list(
  si_distr = gensars, t_start = t_start, t_end = t_end)))

# Extract outputs
tsars = estsars$R$t_end # end of window
Rsars = estsars$R$`Mean(R)`
RsarsCI = matrix(NA, nrow = 2, ncol = length(Rsars))
RsarsCI[1,] = estsars$R$`Quantile.0.025(R)`; RsarsCI[2,] = estsars$R$`Quantile.0.975(R)`

######################################################################
## APE and PMSE solution 
######################################################################
# Outputs take list form [[kbest1, modBest1, kbest2, modBest2]]
# Each modBest is also a list of form [[ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci, alpha, beta, pr]]

# Priors and settings
Rprior = c(1, 5); a = 0.025
# Clean Lam vectors of NAs
Lflu[is.na(Lflu)] = 0; Lsars[is.na(Lsars)] = 0 

# Flu results
Rmodflu = apeEstim(Iflu, genflu, Lflu, Rprior, a, tflu[1], 'flu')
# Best estimates and prediction 
plotAPEWindow(Rmodflu[[2]], 'flu', Rmodflu[[1]], Iplt = Iflu[seq(tflu[1]+1,length(Iflu))], folres)

# Specific 7-day window
Rmodflu7 = apeSpecific(Iflu, genflu, Lflu, Rprior, a, tflu[1], 7)
plotAPEWindow(Rmodflu7[[2]], 'flu7', Rmodflu7[[1]], Iplt = Iflu[seq(tflu[1]+1,length(Iflu))], folres)

# Sars results
Rmodsars = apeEstim(Isars, gensars, Lsars, Rprior, a, tsars[1], 'sars')
# Best estimates and prediction 
plotAPEWindow(Rmodsars[[2]], 'sars', Rmodsars[[1]], Iplt = Isars[seq(tsars[1]+1,length(Isars))], folres)

# Specific 7-day window
Rmodsars7 = apeSpecific(Isars, gensars, Lsars, Rprior, a, tsars[1], 7)
plotAPEWindow(Rmodsars7[[2]], 'sars7', Rmodsars7[[1]], Iplt = Isars[seq(tsars[1]+1,length(Isars))], folres)

######################################################################
## Comparison to EpiEstim
######################################################################

# Compare against EpiEstim outputs for weekly windows as sanity check
Rsars2 = Rmodsars7[[2]][[4]]; Rflu2 = Rmodflu7[[2]][[4]]
Rsars2CI = Rmodsars7[[2]][[5]]; Rflu2CI = Rmodflu7[[2]][[5]]

# Compare APE at weekly windows and EpiEstim
quartz(); par(mfrow=c(2,1))
# Time ranges to plot
ran1 = seq(1, length(tflu)-1); ran2 = seq(1, length(tsars)-1)
# Reprod. num estimates and confidence interval for flu
plot(ran1, Rflu[ran1], type = 'l', bty = 'l', lwd = 2, col='red',
     xlab = paste0("time (k = ", 7, ")"), ylab = 'Rhat flu')
lines(ran1, RfluCI[1, ran1], col = 'red', type = "l", lwd = 1)
lines(ran1, RfluCI[2, ran1], col = 'red', type = "l", lwd = 1)
points(ran1, Rflu2, lwd = 2, col = 'lightgrey')
points(ran1, Rflu2CI[1,], lwd = 2, col = 'lightgrey')
points(ran1, Rflu2CI[2,], lwd = 2, col = 'lightgrey')
# Reprod. num estimates and confidence interval for sars
plot(ran2, Rsars[ran2], type = 'l', bty = 'l', lwd = 2, col='red',
     xlab = paste0("time (k = ", 7, ")"), ylab = 'Rhat sars')
lines(ran2, RsarsCI[1, ran2], col = 'red', type = "l", lwd = 1)
lines(ran2, RsarsCI[2, ran2], col = 'red', type = "l", lwd = 1)
points(ran2, Rsars2, lwd = 2, col = 'lightgrey')
points(ran2, Rsars2CI[1,], lwd = 2, col = 'lightgrey')
points(ran2, Rsars2CI[2,], lwd = 2, col = 'lightgrey')

