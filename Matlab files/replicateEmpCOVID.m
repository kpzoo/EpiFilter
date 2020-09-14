% Renewal estimation and smoothing with recursive filters on empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - all predictions are causal
% - input of empirical data from WHO on COVID-19
% - includes solutions to smoothing and filtering problems
% - no comparison to APE or EpiEstim and naively uses cases
% - assumes serial interval from Ferguson et al

% Directory and if saving
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Extract empirical data from WHO and setup renewal model

% Folder for saving and loading
saveFol = 'Results/COVID/'; loadFol = 'Data/COVID/';
% Possible countries of interest 
countries = {'New Zealand', 'Croatia', 'Greece'};
% Dates of lockdown and release for each
lockdowns = {'26-03-20', '18-03-20', '23-03-20'};
releases = {'14-05-20', '19-04-20', '04-05-20'};

% Decide country to investigate
scenNo = 1; scenNam = countries{scenNo};
disp(['Examining data from ' scenNam]);
% Identifier for saving
namstr = ['_' scenNam];

% Get specific lockdown/release for country
lock = lockdowns{scenNo}; relax = releases{scenNo};
lock = datetime(lock, 'InputFormat', 'dd-MM-yy');
relax = datetime(relax, 'InputFormat', 'dd-MM-yy');

% Load main data and select country
cd(loadFol);

% File with epidemic curves
file = dir('WHO*'); file = file.name;
% Data for all countries
data = readtable(file);

% Select parts with country of interest
id = find(strcmp(data.Country, scenNam));
% Check data
if isempty(id)
    error('Incorrect or unavailable country');
end

% Incidence and death time-series
Iday = data.New_cases(id); nday = length(Iday);
% Dates and deaths
dateRep = data.Date_reported(id); Dday = data.New_deaths(id);
clear('data');

% Reorder dates to be in common format
dateRep = datestr(dateRep); dateRep = datetime(dateRep);
dateRep = datetime(dateRep, 'Format', 'dd-MM-yy'); 
% Get the index of lockdown and release
idlock = find(lock == dateRep); idrelax = find(relax == dateRep);
if isempty(idlock) || isempty(idrelax)
    error('Lockdown and release dates not in data');
end

% Numeric value associated with date, 2000 is pivot
datex = linspace(dateRep(1), dateRep(end), nday);
xval = datenum(datex);
cd(thisDir);

% Gamma serial interval distribution from Ferguson 
distvals.type = 2;
distvals.pm = (1/0.65)^2; distvals.omega = 6.5;
% Serial distribution over all days
serial = serialDistrTypes(nday, distvals);
Pomega = serial(1/distvals.omega);

% Convert to row vectors and times
Iday = Iday'; Dday = Dday'; tday = 1:nday;

% Total infectiousness
Lam = zeros(size(Iday)); 
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = sum(Iday(i-1:-1:1).*Pomegat);    
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
[Rmed, Rlow, Rhigh, Rmean, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, nday, p0, Lam, Iday);

% EpiFilter one-step-ahead predictions 
[predF, predIntF] = recursPredict(Rgrid, pR, Lam, Rmean);

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1 = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(pR(i, :)); prL1(i) = Rcdf(id1);
end

%% Recursive smoother and predictions

% EpiSmoother estimates for single trajectory
[RmedS, RlowS, RhighS, RmeanS, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiSmoother one-step-ahead predictions 
[predS, predIntS] = recursPredict(Rgrid, qR, Lam, RmeanS);

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1S = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(qR(i, :)); prL1S(i) = Rcdf(id1);
end


%% Visualisation and saving


% Bayesian estimates with incidence
figure;
ax(1) = subplot(2, 1, 1);
plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on;
plotCIRaw(tday', Rmed', Rlow', Rhigh', 'r');
plot(tday, Rmean, 'b--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\tilde{R}_s$', 'FontSize', 18)
xlim([tday(1) tday(end)]);
ax(2) = subplot(2, 1, 2);
scatter(tday(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tday(2:end)', predF', predIntF(:, 1), predIntF(:, 2), 'r');
grid off; box off; hold off;
xlim([tday(1) tday(end)]);
ylabel(['$\tilde{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x');
% % Convert from time to dates
% xL = ax(1).XTickLabel; xD = cellfun(@str2num, xL);
% dateVal = date(xD); ax(1).XTickLabel = dateVal;

% Bayesian estimates with incidence
figure;
ax(1) = subplot(2, 1, 1);
plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on;
plotCIRaw(tday', RmedS', RlowS', RhighS', 'r');
plot(tday, RmeanS, 'b--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([tday(1) tday(end)]);
ax(2) = subplot(2, 1, 2);
scatter(tday(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tday(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'r');
grid off; box off; hold off;
xlim([tday(1) tday(end)]);
ylabel(['$\hat{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x');

% Probability at most 1 from EpiFilter
figure
yyaxis left
plot(tday, Iday, '-', 'Color', grey1, 'LineWidth', 2); 
ylabel('$I_s$', 'FontSize', 18)
h = gca; h.YColor = 'k'; 
yyaxis right
hold on;
plot(tday, prL1, 'r-', 'LineWidth', 2);
plot(tday, prL1S, 'b-', 'LineWidth', 2);
h = gca; h.YColor = 'k'; h.YLim(2) = 1.05;
grid off; box off; hold off;
ylabel('$p(R_s \leq 1)$', 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(1) tday(end)]);
if saveTrue
    cd(saveFol);
    saveas(gcf, ['probR1' namstr], 'fig');
    cd(thisDir);
end

% Compare filter and smoother estimates
figure;
ax(1) = subplot(2, 1, 1);
plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on;
plotCIRaw(tday', Rmed', Rlow', Rhigh', 'r');
plot(tday, Rmean, 'b--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\tilde{R}_s$', 'FontSize', 18)
xlim([tday(1) tday(end)]);
ax(2) = subplot(2, 1, 2);
plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on;
plotCIRaw(tday', RmedS', RlowS', RhighS', 'r');
plot(tday, RmeanS, 'b--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([tday(1) tday(end)]);
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'xy');

% Compare filter and smoother predictions
figure;
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predF', predIntF(:, 1), predIntF(:, 2), 'b');
plotCIRaw(xval(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'r');
grid off; box off; hold off;
xlim([xval(2) xval(end)]);
ylabel(['$\tilde{I}_s, \, \hat{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);

%% Publishable figures and saving data

% Bayesian smoothed R against lockdown time
figure;
ax(1) = subplot(2, 1, 1);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on; hdate = gca;
plotCIRaw(xval', RmeanS', RlowS', RhighS', 'r');
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% Predicted incidence
ax(2) = subplot(2, 1, 2);
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'b');
grid off; box off; hold off;
xlim([xval(2) xval(end)]);
ylabel(['$\hat{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x'); datetick('x','dd-mm', 'keeplimits');
% Inset of prob(R <= 1)
axIn = axes('Position',[0.55 0.2 0.3 0.2]);
plot(xval, prL1S, 'r-', 'LineWidth', 2);
hold on;
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
hold off; ylabel('P$(R_s \leq 1)$', 'FontSize', 18)
ylim([0 1.05]); box off; grid off;
datetick('x','dd-mm', 'keeplimits');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['RIsmooth' namstr], 'fig');
    cd(thisDir);
end

% Bayesian filtered R against lockdown time
figure;
ax(1) = subplot(2, 1, 1);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on; hdate = gca;
plotCIRaw(xval', Rmean', Rlow', Rhigh', 'r');
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\tilde{R}_s$', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% Predicted incidence
ax(2) = subplot(2, 1, 2);
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predF', predIntF(:, 1), predIntF(:, 2), 'b');
grid off; box off; hold off;
xlim([xval(2) xval(end)]);
ylabel(['$\tilde{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x'); datetick('x','dd-mm', 'keeplimits');
% Inset of prob(R <= 1)
axIn = axes('Position',[0.55 0.2 0.3 0.2]);
plot(xval, prL1, 'r-', 'LineWidth', 2);
hold on;
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
hold off; ylabel('P$(R_s \leq 1)$', 'FontSize', 18)
ylim([0 1.05]); box off; grid off;
datetick('x','dd-mm', 'keeplimits');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['RIfil' namstr], 'fig');
    cd(thisDir);
end


% Bayesian smoothed and filtered R and filtered predictions
figure;
ax(1) = subplot(2, 1, 1);
hold on; hdate = gca;
plotCIRaw(xval', Rmean', Rlow', Rhigh', [0.94 0.94 0.94]);
plotCIRaw(xval', RmeanS', RlowS', RhighS', 'r');
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$, \, $\tilde{R}_s$', 'FontSize', 18);
ylim([0 5]); xlim([xval(2) xval(end)]); 
datetick('x','dd-mm', 'keeplimits');
% Predicted incidence
ax(2) = subplot(2, 1, 2);
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'b');
grid off; box off; hold off;
xlim([xval(2) xval(end)]);
ylabel(['$\hat{I}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x'); datetick('x','dd-mm', 'keeplimits');
% Inset of prob(R <= 1)
axIn = axes('Position',[0.55 0.2 0.3 0.2]);
plot(xval, prL1S, 'r-', 'LineWidth', 2);
hold on;
plot(xval(idlock)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), hdate.YLim, '--', 'Color', grey2, 'LineWidth', 2);
hold off; ylabel('P$(R_s \leq 1)$', 'FontSize', 18)
ylim([0 1.05]); box off; grid off;
datetick('x','dd-mm', 'keeplimits');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['RIsmoothfil' namstr], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    clear('pRup', 'pstate');
    save(['proc' namstr '.mat']);
    cd(thisDir);
end



