% Given look-back k, run one step ahead Bayesian predictions
function [pred, predInt, prob, R, RInt, pm1, pm2] = getNegBinEmpiricalExcess(k, tPts, I, Lam, trunc)

% Assumptions and notes
% - includes truncation from true start point
% - fixes incorrectly specified first values < k as compared to EpiEstim
% - extra outputs to observe distribution shapes
% - works with empirical distributions from EpiEstim R package
% - posterior for R based on k days look-back (k+1) window length
% - if under k values available just use what is available

% Gamma prior parameters (a = shape, b = scale)
a = 1; b = 2; % from Cori 2013 (1, 5) - used for simulation

% Grouped incidence and infectiousness
B = zeros(1, tPts-1); A = B;

% At each time (tPts) compute historical R, predict at tPts+1
for i = 1:tPts-1
    % Look-back window of k (or less)
    idback = i:-1:max(i-k, 1);
    % Relevant incidence sum (B) and total infectiousness sum (A)
    B(i) = sum(I(idback)); A(i) = sum(Lam(idback));
end

% Range of index time points - starts from trunc
ir = trunc:tPts-1;
% Confidence intervals
predInt = zeros(2, length(ir)); RInt = predInt;

% Parameters of posterior gamma on R
num = a + B;
den = 1/b + A;

% Truncation of related vectors
num = num(ir); den = den(ir);

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

% Distribution parameters to output
pm1 = num; pm2 = 1-p;







