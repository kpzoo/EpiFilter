######################################################################
## Bayesian recursive prediction via EpiFilter
# From: Parag, KV, (2020) “Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves” BioRxiv.
######################################################################

# Notes and assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - can apply causal or smoothing posteriors over R to predict incidence
# - need to set maximum on possible Igrid to interrogate (currently 800)

# Inputs - grid on reproduction numbers (Rgrid), posterior on R (pR), 
# total infectiousness (Lday), mean R esimate using pR (Rmean), confidence level (a and 50%)

# Output - mean prediction (pred) and confidence intervals (predInt)

recursPredictAdj <- function(Rgrid, pR, Lday, Rmean, a, Imax){
  
  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  # Test lengths of inputs
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }

  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]
  
  # Discrete space of possible predictions
  Igrid = 0:Imax; lenI = length(Igrid);
  
  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop("Epidemic size too large")  
  }
  
  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)
  
  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Lday[i]*Rgrid
    # Prob of any I marginalised over Rgrid
    pI = rep(0, lenI)
    
    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # Raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of R
      pI[j] = sum(pIset*pR[i, ])
    }
    
    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)
    
    # Assign prediction results
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]
  }
  # Main outputs: mean and 95% predictions
  recursPredictAdj = list(pred, predInt)
}