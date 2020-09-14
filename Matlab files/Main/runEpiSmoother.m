% Bayesian recursive smoother for renewal models
function [Rmed, Rlow, Rhigh, Rm, qR] = runEpiSmoother(Rgrid, nPts, nday, pR, pRup, pstate)

% Assumptions and notes
% - discrete smoother on grid from Rmin to Rmax
% - assumes runEpiFilter outputs available
% - forward-back algorithm for smoothing applied 

% Last smoothed distribution same as filtered
qR = zeros(nday, nPts); qR(end, :) = pR(end, :);

% Main smoothing equation iteratively computed
for i = nday-1:-1:1
    % Remove zeros
    pRup(i+1, pRup(i+1, :) == 0) = 10^-8;
    
    % Integral term in smoother
    integ = qR(i+1, :)./pRup(i+1, :);
    integ = integ*pstate;
  
    % Smoothed posterior over Rgrid
    qR(i, :) = pR(i, :).*integ;
    % Force a normalisation
    qR(i, :) = qR(i, :)/sum(qR(i, :));
end

% Mean, median and confidence on R
Rm = zeros(1, nday); Rmed = Rm; Rlow = Rm; Rhigh = Rm;
for i = 1:nday
    % Posterior mean and CDF
    Rm(i) = qR(i, :)*Rgrid';
    Rcdf = cumsum(qR(i, :));
    
    % Quantiles from estimates
    Rmed(i) = Rgrid(find(Rcdf > 0.5, 1, 'first'));
    Rlow(i) = Rgrid(find(Rcdf > 0.025, 1, 'first'));
    Rhigh(i) = Rgrid(find(Rcdf > 0.975, 1, 'first'));
end
