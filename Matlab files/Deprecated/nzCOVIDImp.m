% Renewal estimation and smoothing with recursive filters on empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - input of empirical data New Zealand ministry of health
% - solutions to smoothing and filtering problems
% - assumes serial interval from Ferguson/Nishiura et al
% - imported cases distinguished from local ones

% Directory and if saving
thisDir = cd; saveTrue = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Extract empirical data from WHO and setup renewal model

% Folder for saving and loading
saveFol = 'Results/COVID/'; loadFol = 'Data/COVID/';
% Possible countries of interest 
countries = {'New Zealand'};

% Dates of lockdown and release for each
lockdowns = {'26-03-20'}; releases = {'14-05-20'};
% Other dates of interest
secondwave = {'14-08-20'}; present = {'07-10-2020'};

% Decide country to investigate
scenNo = 1; scenNam = countries{scenNo};
disp(['Examining data from ' scenNam]);
% Identifier for saving and data type
namstr = ['_' scenNam]; dataType = 2;

% Get specific lockdown/release for country
lock = lockdowns{scenNo}; relax = releases{scenNo};
lock = datetime(lock, 'InputFormat', 'dd-MM-yy');
relax = datetime(relax, 'InputFormat', 'dd-MM-yy');

% Specific related dates
secwv = secondwave{scenNo}; pres = present{scenNo};
secwv = datetime(secwv, 'InputFormat', 'dd-MM-yy');
pres = datetime(pres, 'InputFormat', 'dd-MM-yy');

% Direct NZ data from government
cd(loadFol);
% File with epidemic curves
file = dir('covid*'); file = file.name;
% Data for New Zealand only
data = readtable(file);
dates = data.DateNotifiedOfPotentialCase;

% Length of time-series
tdates = min(dates):max(dates);
nday = length(tdates); tday = 1:nday;
% Numeric value associated with date, 2000 is pivot
xval = datenum(tdates);

% Collect dates into an epidemic curve
Iday = zeros(1, nday);
for i = 1:nday
    Iday(i) = length(find(dates == tdates(i)));
end

% Possible imported cases
importCase = data.OverseasTravel;
importCase = strcmp(importCase, 'Yes')';
datesImp = dates(importCase); datesLoc = dates(~importCase);

% Collect dates into local vs imported
Iloc = zeros(1, nday); Iimp = Iloc;
for i = 1:nday
    Iloc(i) = length(find(datesLoc == tdates(i)));
    Iimp(i) = length(find(datesImp == tdates(i)));
end
Iloc = Iday; % use to test no imports

cd(thisDir);

% Get the index of lockdown and release
idlock = find(lock == tdates); idrelax = find(relax == tdates);
idsecwv = find(secwv == tdates); idpres = find(pres == tdates);

% Gamma serial interval distribution options <--- need to normalise
distvals.type = 2; SIchoice = 1;
switch(SIchoice)
    case 1
        % From Ferguson et al 2020
        omega = 6.5; pm = (1/0.65)^2;
    case 2
        % From LSHTM
        omega = 3.6; pm = (omega^2)/(3.1^2);
    case 3
        % From Prete
        omega = 2.97; pm = (omega^2)/(3.29^2);
        
end
distvals.pm = pm; distvals.omega = omega;

% Serial distribution over all days
serial = serialDistrTypes(nday, distvals);
Pomega = serial(1/distvals.omega);

% Total infectiousness still based on all cases
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

% EpiFilter estimates using local cases
[Rmed, Rlow, Rhigh, Rmean, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, nday, p0, Lam, Iloc);

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

%% Separe into 4 dates and plot perspectives

% Define indices and portions of epidemic curve
idends = [1 idlock idrelax idsecwv idpres];
nends = length(idends)-1; idset = cell(1, nends); 
% Store relevant variables
Iset = idset; Lset = idset; tset = idset; 
Rmedset = idset; Rlowset = idset; Rhighset = idset;
predset = idset; predIntset = idset; xset = idset;

% For each portion compute estimates and predictions
for i = 1:nends
    % Portion considered 
    idset{i} = 1:idends(i+1);
    Iset{i} = Iloc(idset{i}); Lset{i} = Lam(idset{i});
    % Time series here
    xset{i} = xval(idset{i}); tset{i} = tdates(idset{i});
    
    % Run EpiFilter
    [~, ~, ~, ~, pRset, pRupset, pstateset] = runEpiFilterSm(Rgrid, m, eta,...
        length(Iset{i}), p0, Lset{i}, Iset{i});
    [Rmedset{i}, Rlowset{i}, Rhighset{i}, Rmeanset, qRset] = runEpiSmoother(Rgrid,...
        m, length(Iset{i}), pRset, pRupset, pstateset);

    % One-step-ahead predictions
    [predset{i}, predIntset{i}] = recursPredict(Rgrid, qRset,Lset{i}, Rmeanset);
end

% Compare all estimates of R but over same time-frame
figure;
for i = 1:nends
    subplot(nends, 1, i);
    plotCIRaw(xset{i}', Rmedset{i}', Rlowset{i}', Rhighset{i}', 'r');
    grid off; box off; 
    ylabel('$\hat{R}_s$', 'FontSize', 18);
    xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
    if i == nends
        xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 18);
    end
end

% Compare all predictions of I but over same time-frame
figure;
for i = 1:nends
    subplot(nends, 1, i);
    plotCIRaw(xset{i}(2:end)', predset{i}', predIntset{i}(:, 1), predIntset{i}(:, 2), 'b');
    hold on; 
    scatter(xset{i}', Iset{i}, 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    grid off; box off; hold off;
    ylabel('$\hat{I}_s$', 'FontSize', 18);
    xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
    if i == nends
        xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 18);
    end
end

% Compare all estimates/predictions on same plots
figure; cols = {'b', 'r', 'g', 'c'};
subplot(2, 1, 1); i = 1; 
plotCIRaw(xset{i}', Rmedset{i}', Rlowset{i}', Rhighset{i}', cols{i});
hold on;
for i = 2:nends
    plotCIRaw(xset{i}', Rmedset{i}', Rlowset{i}', Rhighset{i}', cols{i});
end
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');

subplot(2, 1, 2); i = 1; 
plotCIRaw(xset{i}(2:end)', predset{i}', predIntset{i}(:, 1), predIntset{i}(:, 2), cols{i});
hold on;
for i = 2:nends
    plotCIRaw(xset{i}(2:end)', predset{i}', predIntset{i}(:, 1), predIntset{i}(:, 2), cols{i});
end
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
grid off; box off; hold off;
ylabel('$\hat{I}_s$', 'FontSize', 18);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 18);

%% Visualisation and saving

% Compare filter and smoother estimates and predictions
figure;
subplot(2, 1, 1);
plotCIRaw(xval', Rmed', Rlow', Rhigh', 'b');
hold on;
plotCIRaw(xval', RmedS', RlowS', RhighS', 'r');
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('$\tilde{R}_s, \, \hat{R}_s$', 'FontSize', 18);
subplot(2, 1, 2);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predF', predIntF(:, 1), predIntF(:, 2), 'b');
plotCIRaw(xval(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'r');
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('$\tilde{I}_s, \, \hat{I}_s$', 'FontSize', 18);
xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 18);

% Bayesian smoothed R against lockdown time
figure;
ax(1) = subplot(2, 1, 1);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 2);
hold on; 
plotCIRaw(xval', RmeanS', RlowS', RhighS', 'r');
plot(xval(idlock)*ones(1, 2), [min(RlowS) max(RhighS)], '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), [min(RlowS) max(RhighS)], '--', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off;
ylabel('$\hat{R}_s$', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% Predicted incidence
ax(2) = subplot(2, 1, 2);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', predS', predIntS(:, 1), predIntS(:, 2), 'b');
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('$\hat{I}_s$', 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);
linkaxes(ax, 'x'); datetick('x','dd-mm', 'keeplimits');
% Inset of prob(R <= 1)
ax(3) = axes('Position',[0.55 0.22 0.3 0.2]);
plot(xval, prL1S, 'r-', 'LineWidth', 2);
hold on;
plot(xval(idlock)*ones(1, 2), [0 1.05], '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), [0 1.05], '--', 'Color', grey2, 'LineWidth', 2);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('P$(R_s \leq 1)$', 'FontSize', 18); ylim([0 1.05]); 
if saveTrue
    cd(saveFol);
    saveas(gcf, ['RIsmooth' namstr], 'fig');
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



