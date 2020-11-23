function [Rest, Ipred, prL1S] = allFilSmooth(Rgrid, m, eta, nday, p0, Lam, Iloc)

% Assumptions and notes
% - runs filter and smoother in one go

% EpiFilter estimates using local cases
[RmedF, RlowF, RhighF, RmeanF, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, nday, p0, Lam, Iloc);

% EpiFilter one-step-ahead predictions 
[predF, predIntF] = recursPredictMax(Rgrid, pR, Lam, RmeanF, max(Iloc));
predlowF = predIntF(:,1)'; predhighF = predIntF(:,2)';

% EpiSmoother estimates for single trajectory
[RmedS, RlowS, RhighS, RmeanS, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiSmoother one-step-ahead predictions 
[predS, predIntS] = recursPredictMax(Rgrid, qR, Lam, RmeanS, max(Iloc));
predlowS = predIntS(:,1)'; predhighS = predIntS(:,2)';

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1S = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(qR(i, :)); prL1S(i) = Rcdf(id1);
end

% Output data structures for R estimates
Rest.med = [RmedF' RmedS'];
Rest.low = [RlowF' RlowS'];
Rest.high = [RhighF' RhighS'];

% Output data structures for I predictions
Ipred.med = [predF' predS'];
Ipred.low = [predlowF' predlowS'];
Ipred.high = [predhighF' predhighS'];