######################################################################
## Compute APE estimate and get confidence 
# From: Parag, KV, and Donnelly, CA. (2019) “Optimising Renewal Models for 
# Real-Time Epidemic Prediction and Estimation” BioRxiv: 835181.
######################################################################

# Assumptions
# - uses Poisson renewal equation
# - gamma prior required on R

# Inputs - incidence curve (Iday), SI distribution (sidistr), total infectiousness (Lday)
# gamma prior R hyperparameters [a b] (Rprior), confidence (a)
# Output - list of best model and best window size for APE and PMSE

apeEstim <- function(Iday, sidistr, Lday, Rprior, a, trunctime, folres){
  
  # Decide to plot
  wantPlot = 0
  
  # No. time points considered
  tday = length(Iday)
  # Possible window lengths
  k = seq(2, tday-trunctime); lenk = length(k)
  print(paste0(c('Windows from:', k[1], 'to', k[lenk]), collapse = ' '))
  
  # APE scores and outpute
  ape = rep(0, lenk); pmse = ape; apeSet = list()
  
  # For every window length compute R estimate and I predictions
  for(i in 1:lenk){
    # Compute APE(k), output is [[ape, pmse, prob, Rhat, Rhatci, Inex, Inexci]]
    apeSet[[i]] = apePredPost(k[i], sidistr, Lday, Iday, Rprior, a, trunctime)
    # Metrics for squared error
    ape[i] = apeSet[[i]][[1]]; pmse[i] = apeSet[[i]][[2]]
  }
  
  # Optimal model and estimates/predictions
  best1 = which.min(ape); kbest1 = k[best1]
  modBest1 = apeSet[[best1]]
  best2 = which.min(pmse); kbest2 = k[best2]
  modBest2 = apeSet[[best2]]
  
  # Plot ape and pmse curve vs k (normalised by maxima)
  if(wantPlot){
    pdf(file=paste0(folres, 'ape_pmse', '.pdf')) 
    
    # Normalisation
    apeN = ape/max(ape); pmseN = pmse/max(pmse)
    # Smallest and largest values
    metMin = min(min(apeN), min(pmseN)); metMax = 1
    
    plot(k, apeN, type = "l", bty = 'l', lwd = 2, col='red',
         xlab = 'window size (k)', ylab = 'metric', ylim = c(metMin, metMax))
    lines(c(kbest1, kbest1), c(metMin,metMax), col = 'red', type = "h", lwd = 2, lty = 'dashed')
    lines(k, pmseN, col = 'blue', lwd = 2)
    lines(c(kbest2, kbest2), c(metMin,metMax), col = 'blue', type = "h", lwd = 2, lty = 'dashed')
    legend('bottomright', legend = c(paste0("APE: k* = ", kbest1), paste0("PMSE: k* = ", kbest2)), lty = c(1,1), col = c('red', 'blue'), lwd = c(2,2))
    dev.off()
  }
  
  # Output is best model and k
  print(paste0(c('Best k [ape pmse] = ', kbest1, 'and', kbest2), collapse = ' '))
  apeEstim = list(kbest1, modBest1, kbest2, modBest2)
  
}