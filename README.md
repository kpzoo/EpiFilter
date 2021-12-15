# EpiFilter
Maximally informed, mean square error optimised estimates of reproduction numbers (R) over time.

Uses Bayesian recursive filtering and smoothing to maximise the information extracted from the incidence data used. 
Takes a forward-backward approach and provides estimates that combine advantages of EpiEstim and the Wallinga-Teunis method.
Method is exact (and optimal given a grid over R) and deterministic (produces the same answer on the same data).

For full details see: 
Parag KV. (2021) Improved estimation of timevarying reproduction numbers at low case incidence and between epidemic waves. PLoS
Comput. Biol. 17, e1009347. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009347

Here we provide Matlab and R code to implement the main methods described in the text.

Main functions: epiFilter (or epiFilterSm) performs forward filtering to generate causal R estimates; epiSmooth performs backward smoothing to generate R estimates that use all possible information and recursPredict provides one-step-ahead predictions.

Notes on usage:

1) Incidence curve needs to start with a non-zero value
2) Currently only uses gamma serial interval distributions but can be generalised
  - by providing a function for directly computing total infectiousness, Lam
  - e.g. see overall_infectivity function in EpiEstim: https://cran.r-project.org/web/packages/EpiEstim/EpiEstim.pdf
3) Fit of the filtered one-step-ahead predictions gives a measure of model adequacy
  - this follows from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007990


New Zealand branch includes extra functionality for computing local (R, Z) framework described in 
Parag KV, Cowling BJ and Donnelly CA. (2021) Deciphering early-warning signals of SARS-CoV-2 elimination and resurgence from limited data at multiple scales. J. R. Soc. Interface. 18, 20210569. http://doi.org/10.1098/rsif.2021.0569

Here the R of basic EpiFilter is generalised to account for the distinction between local and imported cases. Second, a new metric Z, which assesses
the confidence in local elimination at any time by integrating R is also provided. We apply these to uncover how trends in community transmission in
New Zealand (which features a prolonged period of low or no incidence that destabilises many other R estimators) align with policy action-times. 

Key Matlab code for each case study (New Zealand, Hong Kong and Victoria Australia) are provided in separate branches while the main R code (vignetteImports.R) is in the New Zealand branch.

Further theory surrounding Z and probabilities of elimination can be found here: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008478 
