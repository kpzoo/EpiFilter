% Renewal estimation with recursive Bayesian smoothers on empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - only causal filtering predictions used
% - includes option to filter incidence
% - input of empirical data from EpiEstim R package
% - compares to NB posterior predictions and Bayesian APE
% - discrete replicator/Snyder method employed
% - includes a Brownian motion assumption on R
% - includes smoothing recursive Bayesian method

% Directory, if saving 
thisDir = cd; saveTrue = 1;
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Extract empirical data and Cori estimates

% Folder for saving and loading
saveFol = 'Results/Empirical';  loadFol = 'Data/EpiEstim filter';

% Dataset of choice and if pre-filtering
scenNo = 2; useFilt = 1;

% Define dataset scenarios 
scenNam = {'flu', 'sars'};
scen = scenNam{scenNo};
% Saving name
namstr = ['_' scen '_' num2str(useFilt)];

% Identifier that all files should bear
idfil = [scen '.csv']; idfil2 = [scen 'Filt.csv'];

% Load main data and Cori results
cd(loadFol);

% Get incidence curve data
tday = csvread(['t' idfil]); tday = tday';
IdayNoFilt0 = csvread(['I' idfil]); IdayNoFilt0 = IdayNoFilt0';
% Discretised serial distribution 
wday = csvread(['gen' idfil]); wday = wday';
% Total infectiousness 
LamNoFilt0 = csvread(['L' idfil], 1); 
LamNoFilt0 = [0; LamNoFilt0]; % first val is NA
LamNoFilt0 = LamNoFilt0';

% Cori 7 day estimates (on unfiltered data)
if ~useFilt
    Rcori = csvread(['R' idfil]); Rcori = Rcori';
    RcoriCI = csvread(['RCI' idfil]);
end

% Incidence on moving average
mfil = csvread('mfil.csv');
IdayFilt0 = csvread(['I' idfil2]); IdayFilt0 = IdayFilt0';
LamFilt0 = csvread(['L' idfil2], 1); 
LamFilt0 = [0; LamFilt0]; LamFilt0 = LamFilt0';
cd(thisDir);

% Full time scale is based on Iday
nday0 = length(IdayNoFilt0); tday0 = 1:nday0;

% Truncated time scale is tday
nday = length(tday);
% Truncate Iday and Lam to tday (acts as indices)
Iday1 = IdayNoFilt0(tday); Lam1 = LamNoFilt0(tday);
% Truncate smoothed incidence and Lam
Iday2 = IdayFilt0(tday); Lam2 = LamFilt0(tday);

% Moving average on incidence
if useFilt
    % Filtered incidence
    Iday = Iday2; Lam = Lam2;
    Iday0 = IdayFilt0; Lam0 = LamFilt0;
    disp(['Using filtered incidence with m =' num2str(mfil)]);
else
    % Use raw incidence
    Iday = Iday1; Lam = Lam1;
    Iday0 = IdayNoFilt0; Lam0 = LamNoFilt0;
    disp('No filtering on data');
end

% Space of look-back windows
ks = 1:ceil(nday/2);
nks = length(ks);
disp(['k varies from ' num2str(ks(1)) ' to ' num2str(ks(end))]);


%% APE model selection and prediction

% Posterior incidence predictions
pred = cell(1, nks); predInt = pred;
% Posterior estimates of R over ks
R = pred; RInt = pred;
% APE metric and PMSE
prob = pred; ape = zeros(1, nks); pmse = ape;
% Parameters of NB predictive distribution
pm1 = pred; pm2 = pred;

for i = 1:nks
    % One-step-ahead Bayesian posterior prediction for k
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, pm1{i}, pm2{i}] = ...
        getNegBinEmpiricalExcess(ks(i), nday0, Iday0, Lam0, tday(1));
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
    
    %disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
end

% Best models according to metrics
[apeMin, apeMod] = min(ape);
[pmseMin, pmseMod] = min(pmse);

% Best ks and nGrps
modID = [apeMod pmseMod];
kbest = ks(modID);

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


%% Recursive filter and predictions

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);

% Uniform prior over grid of size m
m = 2000; p0 = (1/m)*ones(1, m);
% Delimited grid defining space of R
Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter estimates for single trajectory
[RmedF, RlowF, RhighF, RmeanF, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, nday, p0, Lam, Iday);

% EpiSmoother estimates for single trajectory
[Rmed, Rlow, Rhigh, Rmean, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiFilter one-step-ahead predictions 
[predF, predIntF] = recursPredict(Rgrid, pR, Lam, RmeanF);
% Smoothed one-step-ahead predictions
[predS, predIntS] = recursPredict(Rgrid, qR, Lam, Rmean);

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1 = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(qR(i, :)); prL1(i) = Rcdf(id1);
end

%% Incidence predictions, R estimates, main figures

% For plotting set k to window length = k + 1
ks = ks + 1; kbest = kbest + 1;
disp(['Bayes k: [ape pmse] = [' num2str(kbest) ']' ]);

% Optimal R window and weekly window
Ropt = R{modID(1)}; id7 = find(ks == 7); R7 = R{id7};
% Plot from second time onwards
tplt = tday(2:end); Iplt = Iday(2:end);

% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
plotCIRaw(tplt', Ropt', RInt{apeMod}(1, :)', RInt{apeMod}(2, :)', 'r');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
subplot(2, 1, 2);
plotCIRaw(tplt', pred{apeMod}', predInt{apeMod}(1, :)', predInt{apeMod}(2, :)', 'r');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
ylabel(['$k^* = $ ' num2str(kbest(1))]);


% APE against predictive error
figure;
yyaxis left
plot(ks, ape, '.', 'linewidth', 2, 'markersize', 20);
hold on; h = gca;
plot([kbest(1) kbest(1)], h.YLim, 'k--', 'LineWidth', 2);
plot([7 7], h.YLim, '--', 'Color', grey2, 'LineWidth', 2);
hold off; 
h = gca; h.YColor = h.XColor;
ylabel('APE', 'FontSize', 18);
yyaxis right
plot(ks, percMiss, '.', 'linewidth', 2, 'markersize', 20);
h = gca; h.YColor = h.XColor;
grid off; box off;
ylabel('$\%$ missed', 'FontSize', 18);
xlabel('$k$ (days)', 'FontSize', 18);
xlim([ks(1)-1 ks(end)+1]);

% Bayesian estimates with incidence
figure;
subplot(2, 1, 1);
hold on;
plotCIRaw(tday', Rmed', Rlow', Rhigh', 'b');
plot(tday, Rmean, 'r--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([tday(2) tday(end)]);
subplot(2, 1, 2);
plotCIRaw(tplt', predF', predIntF(:, 1), predIntF(:, 2), 'r');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
ylabel('$\tilde{I}_s$', 'FontSize', 18)

% Probability at most 1
figure
yyaxis left
plot(tday, Iday, 'b', 'LineWidth', 2); 
ylabel('$I_s$', 'FontSize', 18)
h = gca; h.YColor = 'k'; 
yyaxis right
plot(tday, prL1, 'r', 'LineWidth', 2);
h = gca; h.YColor = 'k'; h.YLim(2) = 1.05;
grid off; box off; 
ylabel('$p(R_s \leq 1)$', 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(1) tday(end)]);

% Compare APE and Bayesian estimates
figure; clear ax
ax(1) = subplot(2, 1, 1);
plot(tplt, ones(size(tplt)), '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tplt', Ropt', RInt{apeMod}(1, :)', RInt{apeMod}(2, :)', 'r');
grid off; box off; hold off;
ylabel(['$\tilde{R}_{\tau(s)} | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(2, 1, 2);
plot(tplt, ones(size(tplt)), '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', RmeanF', RlowF', RhighF', 'c');
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
grid off; box off; hold off;
ylabel('$\tilde{R}_s$ and $\hat{R}_s$', 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
linkaxes(ax, 'xy');

% Compare smoothed vs filtered
figure;
subplot(2, 1, 1);
plot(tplt, ones(size(tplt)), '--', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', RmeanF', RlowF', RhighF', 'c');
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
grid off; box off; hold off;
ylabel('$\tilde{R}_s$ and $\hat{R}_s$', 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tplt(1) tplt(end)]);

subplot(2, 1, 2);
plotCIRaw(tplt', predF', predIntF(:, 1), predIntF(:, 2), 'c');
hold on;
plotCIRaw(tplt', predS', predIntS(:, 1), predIntS(:, 2), 'r');
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel('$\tilde{I}_s$ and $\hat{I}_s$', 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['smoothFilCheck_' namstr], 'fig');
    cd(thisDir);
end

%% Publishable figure and data

% Optimal prediction against 7 day one and smoothing
figure;
% R estimates
ax(1) = subplot(3, 2, 1);
plotCIRaw(tplt', Ropt', RInt{apeMod}(1, :)', RInt{apeMod}(2, :)', 'r');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$\tilde{R}_{\tau(s)} | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(3, 2, 3);
plotCIRaw(tplt', R7', RInt{id7}(1, :)', RInt{id7}(2, :)', 'r');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$\tilde{R}_{\tau(s)} | k = $ ' num2str(ks(id7))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
ax(3) = subplot(3, 2, 5);
hold on;
%plotCIRaw(tday', RmeanF', RlowF', RhighF', 'c');
plotCIRaw(tday', Rmean', Rlow', Rhigh', 'r');
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$\hat{R}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
linkaxes(ax, 'xy');
% I predictions
ax(1) = subplot(3, 2, 2);
plotCIRaw(tplt', pred{apeMod}', predInt{apeMod}(1, :)', predInt{apeMod}(2, :)', 'b');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel(['$\tilde{I}_{\tau(s)} | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(3, 2, 4);
plotCIRaw(tplt', pred{id7}', predInt{id7}(1, :)', predInt{id7}(2, :)', 'b');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel(['$\tilde{I}_{\tau(s)} | k = $ ' num2str(ks(id7))], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
ax(3) = subplot(3, 2, 6);
plotCIRaw(tplt', predF', predIntF(:, 1), predIntF(:, 2), 'b');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1,...
    'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel(['$\tilde{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['compAll_' namstr], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    clear('pstate', 'pRup');
    cd(saveFol);
    save(['data_' namstr '.mat']);
    cd(thisDir);
end



