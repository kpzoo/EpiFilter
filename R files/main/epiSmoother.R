######################################################################
## Bayesian recursive smoothing via EpiFilter
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - must have run epiFilter first to obtain forward distribution pR
# - method makes a backward pass to generate qR

# Inputs - grid on reproduction numbers (Rgrid), size of grid (m), filtered posterior (pR),
# update pre-filter (pRup), max time (nday), state transition matrix (pstate), confidence (a)

# Output - mean (Rmean), median (Rmed), lower (Rlow) amd upper (Rhigh) quantiles of estimates,
# smoothed posterior over R (qR) which is backwards and forwards

epiSmoother <- function(Rgrid, m, pR, pRup, nday, pstate, a){
  
  # Last smoothed distribution same as filtered
  qR = matrix(0, nday, m); qR[nday, ] = pR[nday, ]
  
  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # Remove zeros
    pRup[i+1, pRup[i+1, ] == 0] = 10^-8
    
    # Integral term in smoother
    integ = qR[i+1, ]/pRup[i+1, ]
    integ = integ%*%pstate
    
    # Smoothed posterior over Rgrid
    qR[i, ] = pR[i, ]*integ
    # Force a normalisation
    qR[i, ] = qR[i, ]/sum(qR[i, ]);
  }
  
  # Mean, median estimats of R
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)
  
  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean and CDF
    Rmean[i] = qR[i, ]%*%Rgrid
    Rcdf = cumsum(qR[i, ])
    
    # Quantiles for estimates
    idm = which(Rcdf >= 0.5); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  epiSmoother = list(Rmed, Rhat, Rmean, qR)
}