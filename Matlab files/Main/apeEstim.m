% Run APE estimation method
function [ks, kbest, Ropt, Iopt, Rconf, Iconf, ape] = apeEstim(nday, Iday, Lam)

% Assumptions and notes
% - APE method for optimising window size
% - expects input of incidence curve from renewal model

% Range of windows for APE
ks = 1:ceil(nday/2); nks = length(ks); 
disp(['Window sizes from ' num2str(ks(1)+1) ' to ' num2str(ks(end)+1)]);

% Posterior incidence predictions
pred = cell(1, nks); predInt = pred;
% Posterior estimates of R over ks
R = pred; RInt = pred; prob = pred;
% APE metric and PMSE
ape = zeros(1, nks); pmse = ape; 

for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, ~] = getNegBinAPE2(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
end

% Check for inadmissible values
if any(isnan(ape)) || any(isinf(ape))
    assignin('base', 'Perr', prob);
    error('APE score inadmissible');
end

% Best models according to metrics
[~, apeMod] = min(ape); [~, pmseMod] = min(pmse);
% Best window lengths from metrics
kbest = ks([apeMod pmseMod]);
disp(['Optimal k: [ape pmse] = [' num2str(kbest) ']' ]);

% Best APE model and prediction
Ropt = R{apeMod}; Iopt = pred{apeMod};
Rconf = RInt{apeMod}; Iconf = predInt{apeMod};