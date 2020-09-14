% APE for prediction and model selection in renewal models
clearvars; clc;
close all; tic;
thisDir = cd;

% Assumptions and notes
% - includes option to filter incidence
% - SARS and Flu considered, Cori picked 7 day windows
% - input of empirical data from EpiEstim R package
% - uses NB posterior predictions and Bayesian APE

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save data and test figs
saveTrue = 0; testFig = 0;
% Folder for saving and loading
saveFol = 'Data/empirical filter'; 
loadFol = 'Data/EpiEstim filter';
% Use filtered incidence data
useFilt = 0;

% Dataset of choice
scenNo = 1;

% Define dataset scenarios 
scenNam = {'flu', 'sars', 'meas', 'pox'};
scen = scenNam{scenNo};
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

% Cori 7 day estimates
Rcori = csvread(['R' idfil]); Rcori = Rcori';
RcoriCI = csvread(['RCI' idfil]);

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
    %[pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, pm1{i}, pm2{i}] = getNegBinEmpirical(ks(i), nday, Iday, Lam);
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, pm1{i}, pm2{i}] = ...
        getNegBinEmpiricalExcess(ks(i), nday0, Iday0, Lam0, tday(1));
    
    % APE and predictive MSE for Bayesian
    ape(i) = -sum(log(prob{i}));
    pmse(i) = mean((Iday(2:end) - pred{i}).^2);
    
    disp(['Completed ' num2str(i) ' of ' num2str(nks)]);
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

% Prob and its sum-log
probS = cell2mat(prob');
logS = log(1./probS);
logS = cumsum(logS');
% Successive APE with time
[~, apeIDS] = min(logS, [], 2);

% Mean prediction error with time
predS = cell2mat(pred');
eS2 = (predS - Iday(2:end)).^2;
% Sum of squares then mean
eS2 = cumsum(eS2');
for i = 1:nday-1
    eS2(i, :) = eS2(i, :)/i;
end
[~, pmseIDS] = min(eS2, [], 2);

% [means samps] = bootstrp(1000, @mean, Iday)

%% Incidence predictions, R estimates, main figures

% For plotting set k to window length = k + 1
ks = ks + 1; kbest = kbest + 1;
disp(['Bayes k: [ape pmse] = [' num2str(kbest) ']' ]);

% Successive window sizes to APE and PMSE
kapeIDS = ks(apeIDS); kpmseIDS = ks(pmseIDS);

% Optimal R window and weekly window
Ropt = R{modID(1)}; Rpmse = R{modID(2)};
id7 = find(ks == 7); R7 = R{id7};

% APE estimate and prediction confidence intervals
e1 = Ropt - RInt{apeMod}(1, :); e1 = e1';
e2 = RInt{apeMod}(2, :) - Ropt; e2 = e2';
e1p = pred{apeMod} - predInt{apeMod}(1, :); e1p = e1p';
e2p = predInt{apeMod}(2, :) - pred{apeMod}; e2p = e2p';

% Cori confidence intervals
e1c = Rcori - RcoriCI(1, :); e1c = e1c';
e2c = RcoriCI(2, :) - Rcori; e2c = e2c';
% Truncate Cori estimates <---- check why <-----------------------------
e1c = e1c(1:end-1); e2c = e2c(1:end-1); Rcori = Rcori(1:end-1);

% Seven day window estimation and prediction intervals
e17 = R7 - RInt{id7}(1, :); e17 = e17';
e27 = RInt{id7}(2, :) - R7; e27 = e27';
e17p = pred{id7} - predInt{id7}(1, :); e17p = e17p';
e27p = predInt{id7}(2, :) - pred{id7}; e27p = e27p';

% Plot from second time onwards
tplt = tday(2:end); Iplt = Iday(2:end);
% Best NB predictive parameters
pmbest1 = pm1{apeMod}; pmbest2 = pm2{apeMod};

if testFig

    % Model selection with APE and PMSE
    figure;
    yyaxis left
    plot(ks, ape,'linewidth', 2);
    hold on;
    plot(kbest(1), apeMin, 'o', 'MarkerSize', 10);
    hold off;
    ylabel('APE');
    yyaxis right
    plot(ks, pmse, 'linewidth', 2);
    hold on;
    plot(kbest(2), pmseMin, 'o', 'MarkerSize', 10);
    hold off;
    ylabel('PMSE');
    grid off; box off;
    xlabel('$k$ (days)');
    
    % PMSE against predictive error
    figure;
    yyaxis left
    plot(ks, pmse, '.-', 'linewidth', 2, 'markersize', 20);
    h = gca; h.YColor = h.XColor;
    yyaxis right
    plot(ks, percMiss, '.-', 'linewidth', 2, 'markersize', 20);
    h = gca; h.YColor = h.XColor;
    grid off; box off;
    legend('PMSE', 'predictive error', 'location', 'best');
    xlabel('$k$ (days)');
    
    % One step ahead predictions and metric (APE vs PMSE)
    figure;
    subplot(2, 2, 1);
    plotCI(tplt', Ropt', e1, e2, 'c');
    grid off; box off; 
    ylabel(['$k^* = $ ' num2str(kbest(1))]);
    xlim([tplt(1) tplt(end)]);
    ylim([0 4]);
    subplot(2, 2, 3);
    plotCI(tplt', pred{apeMod}', e1p, e2p, 'c');
    hold on;
    scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    grid off; box off; hold off;
    xlim([tplt(1) tplt(end)]);
    xlabel('$s$ (days)');
    ylabel(['$k^* = $ ' num2str(kbest(1))]);
    
    subplot(2, 2, 2);
    e1 = Rpmse - RInt{pmseMod}(1, :); e1se = e1se';
    e2 = RInt{pmseMod}(2, :) - Rpmse; e2se = e2se';
    plotCI(tplt', Rpmse', e1se, e2se, 'c');
    grid off; box off;
    ylabel(['$k^* = $ ' num2str(kbest(2))]);
    xlim([tplt(1) tplt(end)]);
    ylim([0 4]);
    subplot(2, 2, 4);
    e1 = pred{pmseMod} - predInt{pmseMod}(1, :); e1sep = e1sep';
    e2 = predInt{pmseMod}(2, :) - pred{pmseMod}; e2sep = e2sep';
    plotCI(tplt', pred{pmseMod}', e1sep, e2sep, 'c');
    hold on;
    scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    grid off; box off; hold off;
    xlim([tplt(1) tplt(end)]);
    xlabel('$s$ (days)');
    ylabel(['$k^* = $ ' num2str(kbest(2))]);
    if saveTrue
        cd(saveFol);
        saveas(gcf, ['PMSEAPE_' scen '_' num2str(nks)], 'fig');
        cd(thisDir);
    end
    
    % Indices of interest
    lenx = 20; x = 1:2*max(Iday);
    idx = unique(round(linspace(1, length(pmbest1), lenx)));
    lenx = length(idx); nb = cell(1, lenx);
    
    % All NB predictive distributions
    figure;
    hold on;
    for i = 1:lenx
        % Centre relative to mean
        [xmean, xvar] = nbinstat(pmbest1(idx(i)), pmbest2(idx(i)));
        % Domain to consider
        %x = 1:xmean+3*sqrt(xvar);
        % NB values
        nb{i} = nbinpdf(x, pmbest1(idx(i)), pmbest2(idx(i)));
        plot(x-xmean, nb{i}, 'linewidth', 2);
    end
    grid off; box off; hold off;
    xlabel('$I_{t+1}$');
    ylabel('$P(I_{t+1})$');
    if saveTrue
        cd(saveFol);
        saveas(gcf, ['distr_' scen '_' num2str(nks)], 'fig');
        cd(thisDir);
    end
    
    % Parameters across time
    figure;
    subplot(2, 1, 1);
    plot(tplt, pmbest1, 'linewidth', 2);
    ylabel('$r$');
    grid off; box off; 
    subplot(2, 1, 2);
    plot(tplt, pmbest2, 'linewidth', 2);
    ylabel('$p$');
    grid off; box off;
    xlabel('$s$ (days)');
end

% One step ahead predictions and metric
figure;
subplot(2, 1, 1);
plotCI(tplt', Ropt', e1, e2, 'c');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
%plot(tplt, Rcori, '--', 'color', 'k', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
%ylim([0 4]);
subplot(2, 1, 2);
plotCI(tplt', pred{apeMod}', e1p, e2p, 'c');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
ylabel(['$k^* = $ ' num2str(kbest(1))]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['predRhat_' scen '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% Compare 7 day window estimates from Cori and this implementation
figure;
subplot(2, 1, 1);
plotCI(tplt', R7', e17, e27, 'c');
hold on;
plotCI(tplt', Rcori', e1c, e2c, 'r');
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
subplot(2, 1, 2);
plotCI(tplt', pred{id7}', e17p, e27p, 'c');
hold on;
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');

% Optimal prediction against 7 day one
figure;
ax(1) = subplot(2, 2, 1);
plotCI(tplt', Ropt', e1, e2, 'c');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(2, 2, 3);
plotCI(tplt', R7', e17, e27, 'c');
hold on;
plot(tplt, ones(size(tplt)), 'k--', 'linewidth', 2);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(id7))]);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
linkaxes(ax, 'xy');

mksz = 20;
ax(1) = subplot(2, 2, 2);
plotCI(tplt', pred{apeMod}', e1p, e2p, 'c');
hold on;
%plot(tplt, Iplt, '.', 'color', grey1, 'linewidth', 2, 'markersize', mksz);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel(['$k^* = $ ' num2str(kbest(1))]);
xlim([tplt(1) tplt(end)]);
ax(2) = subplot(2, 2, 4);
plotCI(tplt', pred{id7}', e17p, e27p, 'c');
hold on;
%plot(tplt, Iplt, '.', 'color', grey1, 'linewidth', 2, 'markersize', mksz);
scatter(tplt, Iplt, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel(['$k = $ ' num2str(ks(id7))]);
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['compRI_' scen '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% APE against predictive error
figure;
yyaxis left
plot(ks, ape, '.', 'linewidth', 2, 'markersize', 20);
hold on; h = gca;
plot([kbest(1) kbest(1)], h.YLim, 'k--', 'LineWidth', 2);
plot([7 7], h.YLim, '--', 'Color', grey2, 'LineWidth', 2);
hold off; 
h = gca; h.YColor = h.XColor;
yyaxis right
plot(ks, percMiss, '.', 'linewidth', 2, 'markersize', 20);
h = gca; h.YColor = h.XColor;
grid off; box off;
%legend('APE', 'predictive error', 'location', 'best');
xlabel('$k$ (days)');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['metric_' scen '_' num2str(nks)], 'fig');
    cd(thisDir);
end

% Successive APE results across time
figure;
hax(1) = subplot(2, 1, 1);
stairs(tplt, kapeIDS, 'c', 'linewidth', 2);
grid off; box off;
ylabel('$k^*$');
xlim([tplt(1) tplt(end)]);
hax(2) = subplot(2, 1, 2);
yyaxis left
plot(tplt, Ropt, '-.', 'color', grey2, 'linewidth', 2);
h = gca; h.YColor = h.XColor;
ylabel('$R_s$'); 
ylim([max(0, min(Ropt)-0.3), max(Ropt)+0.3]);
yyaxis right
scatter(tplt, Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
h = gca; h.YColor = h.XColor;
ylabel('$I_{s+1}$');
grid off; box off;
xlim([tplt(1) tplt(end)]);
xlabel('$s$ (days)');
linkaxes(hax, 'x');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['empSucc_' scen '_' num2str(nks)], 'fig');
    cd(thisDir);
end


% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save([scen '_' num2str(nks) '.mat']);
    cd(thisDir);
end



