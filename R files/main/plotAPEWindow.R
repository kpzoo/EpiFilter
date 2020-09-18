######################################################################
## Plot APE estimate and prediction
# From: Parag, KV, and Donnelly, CA. (2019) “Optimising Renewal Models for 
# Real-Time Epidemic Prediction and Estimation” BioRxiv: 835181.
######################################################################

# Assumptions
# - uses model output (Rmod) from apeEstim.R as [[ape, pmse, prob, Rhat, Rhatci, Inexhat, Inexci]]

# Inputs -  apeEstim.R list i.e. Rmod[[2]] for APE or Rmod[[4]] for PMSE, plot name string,
# best window length so Rmod[[1]] for APE and Rmod[[3]] for PMSE, true incidence Iplt
# Output - .eps plots of best R estimate and one-step-ahead incidence predictions

plotAPEWindow <- function(Rmod, plotname, kbest, Iplt, folres){
  # Extract best I(t+1) and R(t) estimates/predictions
  Rhat = Rmod[[4]]; Rhatci = Rmod[[5]]
  Inexhat = Rmod[[6]]; Inexhatci = Rmod[[7]]
  
  # Check lengths
  if (length(Rhat) != length(Inexhat)){
    print('Inconsistent incidence and reprod. num vectors')
  }else{
    # Length of time (relative)
    tset = 1:length(Rhat)
    
    # Two panel plot of estimates and predictions
    pdf(file=paste0(folres, plotname, '.pdf')) 
    par(mfrow=c(2,1))
    # Reprod. num estimates and confidence interval
    plot(tset, Rhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (k = ", kbest, ")"), ylab = 'reprod. number')
    polygon(c(tset, rev(tset)), c(Rhatci[1,], rev(Rhatci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rhatci[3,], rev(Rhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    
    # Incidence predictions and confidence interval
    plot(tset, Inexhat, type = 'l', bty = 'l', lwd = 2, col='blue',
         xlab = paste0("time (k = ", kbest, ")"), ylab = 'incidence',  ylim = c(0, max(Iplt)+30))
    polygon(c(tset, rev(tset)), c(Inexhatci[1,], rev(Inexhatci[2,])),
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Inexhatci[3,], rev(Inexhatci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    points(tset, Iplt, pch = 19, col = 'gray')
    dev.off()
  
  }
}