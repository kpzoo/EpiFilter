% APE for prediction and model selection in renewal models
function [pMiss, modID, ape, Rmse, Imse] = apeForecastFn(scenNo, ks, nks, tday0, nday0)

% Assumptions and notes
% - all time vectors truncated for plotting
% - uses NB posterior predictions and Bayesian APE
% - APE compares sequential predictions with true values
% - simulates a single epidemic, predicts I(t), estimates R(t)

% Simulate epidemic scenarios and truncate
Iwarn = 1; % ensure no warnings
while Iwarn
    [Iday, Lam, Rtrue, tday, Iwarn] = epiSimAPE(scenNo, tday0, nday0, 1);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period
nday = length(tday);

% Posterior incidence predictions 
pred = cell(1, nks); predInt = pred; 
% Posterior estimates of R over ks
R = pred; RInt = pred;
% APE metric and PMSE
prob = pred; ape = zeros(1, nks); 
for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}] = getNegBinAPE(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));    
end

% Best models according to metrics
[~, modID] = min(ape);
% Best ks and nGrps
kbest = ks(modID); 
disp(['Bayes k: ape = ' num2str(kbest)]);

% Accuracy of predictions and estimates over ks
Rmse = zeros(1, nks); Imse = Rmse; pMiss = Rmse;
% Truncate to actual estimated/predicted times
Itrue = Iday(2:end); Rtrue = Rtrue(2:end); 
for i = 1:nks
    predI = predInt{i};
    % Ids at which incidence outside credible interval
    idout = union(find(Itrue < predI(1, :)), find(Itrue > predI(2, :)));
    % Percentage missed 
    pMiss(i) = 100*length(idout)/length(Itrue);
    % MSE in predictions
    nonz = find(Itrue);
    Imse(i) = mean((1 - pred{i}(nonz)./Itrue(nonz)).^2);
    % MSE in estimation
    Rmse(i) = mean((1 - R{i}./Rtrue).^2);
end


