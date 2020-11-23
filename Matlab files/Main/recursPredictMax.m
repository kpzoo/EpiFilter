% Given discrete R distribution get Poiss predictions
function [pred, predInt] = recursPredictMax(Rgrid, pR, Lam, Rmean, Imax)

% Assumptions and notes
% - adapts its maxima on the incidence space
% - uses posterior over R from recursive filter
% - computes APE score for predictions

% No. points and days and range
[nday, m] = size(pR); ir = 1:nday-1;

% Test lengths of inputs
if length(Rgrid) ~= m || length(Lam) ~= nday
    error('Input vectors of incorrect dimension');
end

% Discrete space of possible predictions
Igrid = 0:round(2*Imax); lenI = length(Igrid);
% Check if close to upper bound
pred0 = Lam.*Rmean; pred0 = pred0(2:end);
if any(pred0 > 0.9*max(Igrid))
    assignin('base', 'predErr', pred0);
    error('Epidemic size too large');
end

% Prediction cdf and quantiles
Fpred = zeros(nday-1, lenI);
predInt = zeros(nday-1, 2);
pred = zeros(1, nday-1);

% At every time construct CDF of predictions
for i = ir
   % Compute rate from Poisson renewal
    rate = Lam(i)*Rgrid;
    
    % Prob of any I marginalised over Rgrid
    pI = zeros(1, lenI);
    
    % Probs of observations 1 day ahead
    for j = 1:lenI
        % Raw probabilities of Igrid
        pIset = poisspdf(Igrid(j), rate);
        % Normalised by probs of R
        pI(j) = sum(pIset.*pR(i, :));
    end
    
    % Quantile predictions and CDF at i+1
    Fpred(i, :) = cumsum(pI)/sum(pI); % added /sum
    idlow = find(Fpred(i, :) >= 0.025, 1, 'first');
    idhigh = find(Fpred(i, :) >= 0.975, 1, 'first');
    
    % Assign prediction results
    predInt(i, 1) = Igrid(idlow); predInt(i, 2) = Igrid(idhigh);
    % Mean prediction by integrating over grid
    if i > 1
        pred(i) = Igrid*(pI/sum(pI))';
    else
        % Initialise with renewal mean 
        pred(i) = pred0(1);
    end
end

