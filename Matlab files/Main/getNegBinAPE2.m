% Given look-back k, run one step ahead Bayesian predictions
function [pred, predInt, prob, R, RInt, wins] = getNegBinAPE2(k, tPts, I, Lam)

% Assumptions and notes
% - vectorised for speed
% - posterior for R based on k days look-back
% - if under k values available just use what is available
% - point and interval prediction from negative binomial 
% - conjugate gamma priors used for renewal model

% Gamma prior parameters (a = shape, b = scale)
a = 1; b = 2; % from Cori 2013 b = 5
% Range of index time points
ir = 1:tPts-1; lenir = length(ir);

% Grouped incidence and infectiousness
B = zeros(1, lenir); A = B;
% Confidence intervals
predInt = zeros(2, lenir); RInt = predInt;

% At each time (tPts) compute historical R, predict at tPts+1
for i = ir
    % Look-back window of k (or less)
    idback = i:-1:max(i-k, 1); 
    % Relevant incidence sum (B) and total infectiousness sum (A)
    B(i) = sum(I(idback)); A(i) = sum(Lam(idback));
end

% Take length of last window (likely to be complete)
wins = length(idback);

% Parameters of posterior gamma on R
num = a + B;
den = 1/b + A;

% Posterior mean of R and its 99% confidence
R = num./den;
RInt(1, :) = gaminv(0.025, num, 1./den);
RInt(2, :) = gaminv(0.975, num, 1./den);

% Predictive point for R (mean)
pred = Lam(ir).*R; % Lam(i) = Lam_{t+1}

% Negative binomial confidence (interval estimate)
p = Lam(ir)./den; p = p./(p + 1);
predInt(1, :) = nbininv(0.025, num, 1-p);
predInt(2, :) = nbininv(0.975, num, 1-p);

% Prob of next incidence value - Pilatowska
prob = nbinpdf(I(ir+1), num, 1-p);






