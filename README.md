# EpiFilter
Maximally informed, mean square error optimised estimates of reproduction numbers (R) over time.

Uses Bayesian recursive filtering and smoothing to maximise the information extracted from the incidence data used. 
Takes a forward-backward approach and provides estimates that combine advantages of EpiEstim and the Wallinga-Teunis method.
Method is exact (and optimal given a grid over R) and deterministic (produces the same answer on the same data).

For full details see: 
Parag KV (2020) “Improved real-time estimation of reproduction numbers at low case incidence and between epidemic waves” medRxiv 2020.09.14.20194589.

Here we provide Matlab and R code to implement the main methods described in the text.

Main functions: epiFilter (or epiFilterSm) performs forward filtering to generate causal R estimates; epiSmooth performs backward smoothing to generate R estimates that use all possible information and recursPredict provides one-step-ahead predictions.

Notes on usage:
1) Incidence curve needs to start with a non-zero value
2) Currently only uses gamma serial interval distributions but can be generalised
3) Fit of the filtered one-step-ahead predictions gives a measure of model adequacy

An application of this approach to New Zealand together with a computation of the confidence in elimination can be found on the New Zealand branch. This extension
separates cases based on whether they are local to the region of interest or imported. 

Details are here:
Parag KV and Donnelly CA (2020) "Novel early-warning system for evaluating control and elimination of SARS-CoV-2 reveals alignment of policy with 
local transmission in New Zealand" medRxiv 2020.11.23.20236968.
