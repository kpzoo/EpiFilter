% Renewal estimation with recursive Bayesian filters
clearvars; clc; close all; tic;

% Assumptions and notes
% - consider cases with small incidence or second waves
% - examine long windows, APEestim and EpiFilter
% - cleaned up version of filterDyingEpi

% Folders for saving data
thisDir = cd; saveFol = 'Results/Dying/scen ';
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Setup simulated epidemic

% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Choose scenario and generation time distribution
scenNo = 1; distNo = 2;

% Specific scenario parameters
switch(scenNo)
    % Rs are distinct reprod nums, ts are switch points
    case 1
        % Rapidly controlled epidemic
        Rch = [2 0.5]; tch = 100;
    case 2
        % Square wave (recovered rapid control)
        Rch = [2.5 0.5 2.5]; tch = [70 230];
    case 3
        % Three stage control with sines
        Rch = [4 0.6 2 0.2]; tch = [40 80 150];
    case 4
        % Exponential rise and fall
        tch = 30; Rch = [];
    case 5
        % Two stage control with noise
        Rch = [4 0.8 0.4]; tch = [15 180];
    case 6
        % Second wave dynamics (sines)
        Rch = [1.3 1.2]; tch = 3;
    case 7
        % Long period of low R between transmission and noise
        Rch = [2.5 0.1 1.3]; tch = [50 200];
    case 8
        % Exponential rise and fall and rise
        tch = [40 190]; Rch = [];
end
simVals.Rch = Rch; simVals.tch = tch;

% Simulate epidemic scenarios and truncate
Iwarn = 1; % ensure no warnings
while Iwarn
    [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimScenDie(scenNo, nday0, distNo, simVals);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period
nday = length(tday);

% Saving data and figs
saveTrue = 0; saveFol = join([saveFol num2str(scenNo)]);
namstr = [num2str(nday) '_' num2str(scenNo) '_' num2str(distNo) '_'];

%% Bayesian recursive filter 

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);

% Uniform prior over grid of size m
m = 2000; p0 = (1/m)*ones(1, m);
Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter method
[Rmed, Rlow, Rhigh, Rm, pR] = runEpiFilter(Rmin, Rmax, m, eta, nday, p0, Lam, Iday);

% One step ahead predicitions from EpiFilter 
[predF, predIntF] = recursPredict(Rgrid, pR, Lam, Rm);

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last');
prL1 = zeros(1, nday); 

% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and critical R
    Rcdf = cumsum(pR(i, :));
    % Probability R <= 1
    prL1(i) = Rcdf(id1);
end

%% APE/PMSE estimate as a reference

% Range of windows for APE
ks = 1:ceil(nday/2); nks = length(ks); 
disp(['Window sizes from ' num2str(ks(1)+1) ' to ' num2str(ks(end)+1)]);

% Posterior incidence predictions
pred = cell(1, nks); predInt = pred;
% Posterior estimates of R over ks
R = pred; RInt = pred; prob = pred;
% APE metric and PMSE
ape = zeros(1, nks); pmse = ape; wins = ape;

for i = 1:nks
    % One step ahead Bayesian posterior prediction for k
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, wins(i)] = getNegBinAPE2(ks(i), nday, Iday, Lam);
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
    
    % Print every 10
    if rem(i, 20) == 0
        disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
    end
end

% Check for inadmissible values
if any(isnan(ape)) || any(isinf(ape))
    assignin('base', 'Perr', prob);
    error('APE score inadmissible');
end

% Best models according to metrics
[apeMin, apeMod] = min(ape); [pmseMin, pmseMod] = min(pmse);
% Best window lengths from metrics
kbest = ks([apeMod pmseMod]);
disp(['Bayes k: [ape pmse] = [' num2str(kbest) ']' ]);

% Accuracy of predictions over ks
percMiss = zeros(1, nks);
Itrue = Iday(2:end);
for i = 1:nks
    predI = predInt{i};
    % Ids at which incidence outside credible interval
    idout = union(find(Itrue < predI(1, :)), find(Itrue > predI(2, :)));
    % Percentage missed
    percMiss(i) = 100*length(idout)/length(Itrue);
end

% Best APE model and prediction
Ropt = R{apeMod}; Iopt = pred{apeMod};
Rconf = RInt{apeMod}; Iconf = predInt{apeMod};

% APE with a long window (month)
klong = 30; disp(['Long window size: ' num2str(klong+1)]);
% Corresponding estimates and predictions
RoptL = R{klong}; IoptL = pred{klong};
RconfL = RInt{klong}; IconfL = predInt{klong};

%% Figures of individual estimates and predictions

% Set k to window length = k + 1 
ks = ks + 1; kbest = kbest + 1; klong = klong+1;
% Plot from second time onwards
tplt = tday(2:end); Rplt = Rtrue(2:end); Iplt = Iday(2:end);

% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
plotCIRaw(tplt', Ropt', Rconf(1, :)', Rconf(2, :)', 'r');
grid off; box off; hold off;
ylabel(['$\tilde{R}_s | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
subplot(2, 1, 2);
hold on;
scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
plotCIRaw(tplt', Iopt', Iconf(1, :)', Iconf(2, :)', 'r');
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
ylabel(['$\tilde{I}_s | k^* = $ ' num2str(kbest(1))], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);

% Long window predictions and metric
figure;
subplot(2, 1, 1);
hold on;
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
plotCIRaw(tplt', RoptL', RconfL(1, :)', RconfL(2, :)', 'r');
grid off; box off; hold off;
ylabel(['$\tilde{R}_s | k^* = $ ' num2str(klong)], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
subplot(2, 1, 2);
hold on;
scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
plotCIRaw(tplt', IoptL', IconfL(1, :)', IconfL(2, :)', 'r');
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
ylabel(['$\tilde{I}_s | k^* = $ ' num2str(klong)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);

% Bayesian estimates with incidence
figure;
subplot(2, 1, 1);
plot(tday, Rtrue, 'k--', 'LineWidth', 2);
hold on;
plotCIRaw(tday', Rmed', Rlow', Rhigh', 'r');
plot(tday, Rm, 'b--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([tday(1) tday(end)]);
subplot(2, 1, 2);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tday(2:end)', predF', predIntF(:, 1), predIntF(:, 2), 'r');
grid off; box off; hold off;
xlim([tday(2) tday(end)]);
ylabel(['$\tilde{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);

% APE against predictive error
figure;
yyaxis left
plot(ks, ape, '.', 'linewidth', 2, 'markersize', 20);
hold on; h = gca;
plot([kbest(1) kbest(1)], h.YLim, 'k--', 'LineWidth', 2);
hold off; 
h = gca; h.YColor = h.XColor;
yyaxis right
plot(ks, percMiss, '.', 'linewidth', 2, 'markersize', 20);
h = gca; h.YColor = h.XColor;
grid off; box off;
%legend('APE', 'predictive error', 'location', 'best');
xlabel('$k$ (days)', 'FontSize', 18);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['metric_' namstr], 'fig');
    cd(thisDir);
end

% Probability at most 1 from EpiFilter
figure
yyaxis left
stairs(tday, Iday, 'b', 'LineWidth', 2); 
ylabel('$I_s$', 'FontSize', 18)
h = gca; h.YColor = 'k'; 
yyaxis right
stairs(tday, prL1, 'r', 'LineWidth', 2);
h = gca; h.YColor = 'k'; h.YLim(2) = 1.05;
grid off; box off; 
ylabel('$p(R_s \leq 1)$', 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(1) tday(end)]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['probR1_' namstr], 'fig');
    cd(thisDir);
end

%% Publishable plots

% Comparison of APE, long window and EpiFilter
figure;
% APE estimates
subplot(3, 2, 1);
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tplt', Ropt', Rconf(1, :)', Rconf(2, :)', 'r');
grid off; box off; hold off;
ylabel(['$\tilde{R}_{\tau(s)} | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
% APE predictions
subplot(3, 2, 2);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tplt', Iopt', Iconf(1, :)', Iconf(2, :)', 'r');
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
ylabel(['$\tilde{I}_{\tau(s)} | k^* = $ ' num2str(kbest(1))], 'FontSize', 18)

% Long window estimates
subplot(3, 2, 3);
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tplt', RoptL', RconfL(1, :)', RconfL(2, :)', 'r');
grid off; box off; hold off;
ylabel(['$\tilde{R}_{\tau(s)} | k = $ ' num2str(klong)], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
% Long window predictions
subplot(3, 2, 4);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tplt', IoptL', IconfL(1, :)', IconfL(2, :)', 'r');
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
ylabel(['$\hat{I}_{\tau(s)} | k = $ ' num2str(klong)], 'FontSize', 18)

% EpiFilter estimates
subplot(3, 2, 5);
plot(tplt, Rplt, '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', Rm', Rlow', Rhigh', 'r');
grid off; box off; hold off;
ylabel(['$\hat{R}_{s} | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)', 'FontSize', 18);
% Long window predictions
subplot(3, 2, 6);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tplt', predF', predIntF(:, 1), predIntF(:, 2), 'r');
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
ylabel(['$\tilde{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['compAll_' namstr 'die'], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Remove unneeded variables and save
clear('pstate', 'Rcdf0', 'pI', 'pRup', 'rate', 'prob');
if saveTrue
    cd(saveFol);
    save(['dieScen' namstr  'data' '.mat']);
    cd(thisDir);
end


