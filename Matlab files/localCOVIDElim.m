% Renewal estimation and smoothing with recursive filters on empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - works with NZ or HK data
% - includes updated figures and elimination times
% - only looks at resurgence and compares imported vs non-imported
% - input of empirical data New Zealand ministry of health
% - solutions to smoothing and filtering problems
% - assumes serial interval from Ferguson/Nishiura et al

% Directory and if saving
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Extract empirical data from WHO and setup renewal model

% Possible countries of interest 
countries = {'New Zealand', 'Hong Kong'};
% Decide country to investigate
scenNo = 2; scenNam = countries{scenNo};
disp(['Examining data from ' scenNam]);

% Dates of lockdown and release for each
lockdowns = {'25-03-20', '29-03-20'}; releases = {'14-05-20', '27-05-20'};
% Other dates of interest
secondwave = {'13-08-20', '20-07-20'}; present = {'07-10-2020', '23-11-20'};

% Identifier for saving and data type
namstr = ['_' scenNam]; dataType = 2;
% Folder for saving and loading
switch scenNo
    case 1
        % NZ data and results
        saveFol = 'Results/NZ/'; loadFol = 'Data/NZ/';
    case 2
        % HK data and results
        saveFol = 'Results/HK/'; loadFol = 'Data/HK/';
    otherwise
        error('Country not specified');
end

% Get specific lockdown/release for country
lock = lockdowns{scenNo}; relax = releases{scenNo};
lock = datetime(lock, 'InputFormat', 'dd-MM-yy');
relax = datetime(relax, 'InputFormat', 'dd-MM-yy');

% Specific related dates
secwv = secondwave{scenNo}; pres = present{scenNo};
secwv = datetime(secwv, 'InputFormat', 'dd-MM-yy');
pres = datetime(pres, 'InputFormat', 'dd-MM-yy');

% Elimination dates
elimination = {'08-06-20', '04-07-20'}; 
elim = elimination{scenNo}; elim = datetime(elim, 'InputFormat', 'dd-MM-yy');

% Data from country of interest
cd(loadFol);

switch scenNo
    case 1
        % File with NZ epidemic curves
        file = dir('covid-cases-23oct20_confirmed.xlsx'); file = file.name;
        % Data for New Zealand only
        data = readtable(file);
        dates = data.DateNotifiedOfPotentialCase;
        
        % Possible imported cases
        importCase = data.OverseasTravel;
        importCase = strcmp(importCase, 'Yes')';
        
        % Local and imported
        datesImp = dates(importCase); datesLoc = dates(~importCase);
        % Starting time date of interest
        start = min(dates);
        
    case 2
        % File with HK case dates
        file = dir('hkcase_20201124.csv'); file = file.name;
        % Data for Hong Kong only
        data = readtable(file);
        
        % Take onset date and find NaNs
        dates = data.onset_date; 
        %dates = data.confirm_date;
        missing = isnan(datenum(dates));
        % Replace missing with confirmation dates
        dates(missing) = data.confirm_date(missing);
        
        % Possible imported cases
        reclass = string(data.reclassification);
        importCase = reclass == "Imported";
        % Include linked cases
        maybeImport = reclass == "Epidemiologically linked with Imported";
        importCase = importCase | maybeImport;
        
        % Get indices and corresponding dates
        datesImp = dates(importCase); datesLoc = dates(~importCase);
         % Starting time date of interest
        start = datetime('23-01-20', 'InputFormat', 'dd-MM-yy');
        
end

% Length of time-series
tdates = start:pres;
nday = length(tdates); tday = 1:nday;
% Numeric value associated with date, 2000 is pivot
xval = datenum(tdates);

% Collect dates into an epidemic curve
Iday = zeros(1, nday);
for i = 1:nday
    Iday(i) = length(find(dates == tdates(i)));
end

% Collect dates into local vs imported
Iloc = zeros(1, nday); Iintro = Iloc;
for i = 1:nday
    Iloc(i) = length(find(datesLoc == tdates(i)));
    Iintro(i) = length(find(datesImp == tdates(i)));
end

cd(thisDir);

% Get the index of lockdown and release
idlock = find(lock == tdates); idrelax = find(relax == tdates);
idsecwv = find(secwv == tdates); idpres = find(pres == tdates);
% Get index of elimination 
idelim = find(elim == tdates);

% Gamma serial interval distribution options 
distvals.type = 2; SIchoice = 1;
disp(['Serial interval choice: ' num2str(SIchoice)]);
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
    case 4
        % From Nishiura
        omega = 4.7; pm = (omega^2)/(2.9^2);
        
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

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);

% Uniform prior over grid of size m
m = 2000; p0 = (1/m)*ones(1, m);
% Delimited grid defining space of R
Rgrid = linspace(Rmin, Rmax, m);

%% Recursive filter/smoother, predictions and elimination

% Transmission esimates when distinguish imports
[Rimp, Iimp, p1imp] = allFilSmooth(Rgrid, m, eta, nday, p0, Lam, Iloc);
% Transmission esimates when ignore imports
[Rtot, Itot, p1tot] = allFilSmooth(Rgrid, m, eta, nday, p0, Lam, Iday);

% Times to interrogate for elimination and probs
tindex = 2:length(Iday); lenind = length(tindex);
zfloc = zeros(1, lenind); zsloc = zfloc; zftot = zfloc; zstot = zfloc;
zflocM = zfloc; zslocM = zfloc; zftotM = zfloc; zstotM = zfloc;

% At every time point compute elimination potential
for i = 1:lenind
    % Prob of elimination when distinguish imports
    [zfloc(i), zsloc(i), zflocM(i), zslocM(i)] = getProbEliminFilterDist(Iloc(1:tindex(i)), Iday(1:tindex(i)),...
        Rgrid, m, eta, distvals, p0);
    % Prob of elimination when ignore imports
    [zftot(i), zstot(i), zftotM(i), zstotM(i)] = getProbEliminFilterDist(Iday(1:tindex(i)), Iday(1:tindex(i)),...
        Rgrid, m, eta, distvals, p0);
end

% Find elimination dates
elimLoc = [find(zfloc >= 0.95, 1, 'first') find(zsloc >= 0.95, 1, 'first')];
elimNaive = [find(zftot >= 0.95, 1, 'first') find(zstot >= 0.95, 1, 'first')];
elimLoc = tdates(elimLoc+1); elimNaive = tdates(elimNaive+1);

%% Publishable figures

% All dates as vertical lines
pltdate = kron(xval([idlock idrelax idsecwv idelim]), ones(1, 2));
pltdate = reshape(pltdate, [2 4]);
% Y plotting for Z
Zyrange = [0, 1.05]; Zyrange = reshape(repmat(Zyrange, [1, 4]), [2 4]);

% Y plotting range for R
Ryrange = [min(Rimp.low(:, 2)), max(Rtot.high(:, 2))];
RyrangeF = [min(Rimp.low(:, 1)), max(Rtot.high(:, 1))];
Ryrange = reshape(repmat(Ryrange, [1, 4]), [2 4]);
RyrangeF = reshape(repmat(RyrangeF, [1, 4]), [2 4]);

% Y plotting range for I
Iyrange = [min(Iimp.low(:, 2)), max(Itot.high(:, 2))];
Iyrange = reshape(repmat(Iyrange, [1, 4]), [2 4]);
IyrangeF = [min(Iimp.low(:, 1)), max(Itot.high(:, 1))];
IyrangeF = reshape(repmat(IyrangeF, [1, 4]), [2 4]);


% Bayesian smoothed R against lockdown time
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
ax(1) = subplot(2, 2, 1);
hold on; 
plotCIRaw(xval', Rimp.med(:, 2), Rimp.low(:, 2), Rimp.high(:, 2), 'r');
plot(pltdate, Ryrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (local)', 'FontSize', 18)
ylabel('smoothed $R_s$ (local)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');

ax(3) = subplot(2, 2, 3);
hold on; 
plotCIRaw(xval', Rtot.med(:, 2), Rtot.low(:, 2), Rtot.high(:, 2), 'r');
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
plot(pltdate, Ryrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (total)', 'FontSize', 18)
ylabel('smoothed $R_s$ (naive)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
linkaxes(ax([1 3]), 'y');
xlabel('time, $s$ (days)', 'FontSize', 18);

% Predicted causal incidence
ax(2) = subplot(2, 2, 2);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Iimp.med(:, 2), Iimp.low(:, 2), Iimp.high(:, 2), 'b');
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (local)', 'FontSize', 18);
ylabel('smoothed $I_s$ (local)', 'FontSize', 18);

ax(4) = subplot(2, 2, 4);
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Itot.med(:, 2), Itot.low(:, 2), Itot.high(:, 2), 'b');
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (total)', 'FontSize', 18);
ylabel('smoothed $I_s$ (naive)', 'FontSize', 18);
xlabel('time, $s$ (days)', 'FontSize', 18);
linkaxes(ax([2 4]), 'y');

if saveTrue
    cd(saveFol);
    saveas(gcf, ['estPred' namstr], 'fig');
    cd(thisDir);
end

% Bayesian filtered R against lockdown time
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
ax(1) = subplot(2, 2, 1);
hold on; 
plotCIRaw(xval', Rimp.med(:, 1), Rimp.low(:, 1), Rimp.high(:, 1), 'r');
plot(pltdate, RyrangeF, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (local)', 'FontSize', 18)
ylabel('filtered $R_s$ (local)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');

ax(3) = subplot(2, 2, 3);
hold on; 
plotCIRaw(xval', Rtot.med(:, 1), Rtot.low(:, 1), Rtot.high(:, 1), 'r');
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
plot(pltdate, RyrangeF, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (total)', 'FontSize', 18)
ylabel('filtered $R_s$ (naive)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
linkaxes(ax([1 3]), 'y');
xlabel('time, $s$ (days)', 'FontSize', 18);

% Predicted causal incidence
ax(2) = subplot(2, 2, 2);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Iimp.med(:, 1), Iimp.low(:, 1), Iimp.high(:, 1), 'b');
plot(pltdate, IyrangeF, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (local)', 'FontSize', 18);
ylabel('filtered $I_s$ (local)', 'FontSize', 18);

ax(4) = subplot(2, 2, 4);
scatter(xval(2:end), Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Itot.med(:, 1), Itot.low(:, 1), Itot.high(:, 1), 'b');
plot(pltdate, IyrangeF, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (total)', 'FontSize', 18);
ylabel('filtered $I_s$ (naive)', 'FontSize', 18);
xlabel('time, $s$ (days)', 'FontSize', 18);
linkaxes(ax([2 4]), 'y');

if saveTrue
    cd(saveFol);
    saveas(gcf, ['estFil' namstr], 'fig');
    cd(thisDir);
end

% Prob of elimination across time
figure('Renderer', 'painters', 'Position', [10 10 1100 500]);
subplot(1, 2, 1);
stairs(xval(2:end), 100*zstot, 'Color', 'g', 'LineWidth', 2);
hold on;
stairs(xval(2:end), 100*zsloc, 'Color', [0.07,0.62,1.00], 'LineWidth', 2);
plot(xval(idlock)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idelim)*ones(1, 2), 100*[0 1.05], '-', 'Color', 'r', 'LineWidth', 2);
plot(xval(2:end), 100*0.95*ones(size(xval(2:end))), 'k--', 'LineWidth', 1);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('smoothed $Z_s$', 'FontSize', 18); ylim(100*[0 1.05]); 
xlabel('time, $s$ (days)', 'FontSize', 18);
legend('naive', 'local', 'Location', 'best'); legend boxoff;   
subplot(1, 2, 2);
stairs(xval(2:end), 100*zftot, 'Color', 'g', 'LineWidth', 2);
hold on;
stairs(xval(2:end), 100*zfloc, 'Color', [0.07,0.62,1.00], 'LineWidth', 2);
plot(xval(idlock)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), 100*[0 1.05], '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(idelim)*ones(1, 2), 100*[0 1.05], '-', 'Color', 'r', 'LineWidth', 2);
plot(xval(2:end), 100*0.95*ones(size(xval(2:end))), 'k--', 'LineWidth', 1);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('filtered $Z_s$', 'FontSize', 18); ylim(100*[0 1.05]); 
xlabel('time, $s$ (days)', 'FontSize', 18);
legend('local', 'naive', 'Location', 'best'); legend boxoff;  

if saveTrue
    cd(saveFol);
    saveas(gcf, ['elim' namstr], 'fig');
    cd(thisDir);
end

% Examine local and total incidence curve
figure;
stairs(xval, Iday, 'b-', 'LineWidth', 2);
hold on;
stairs(xval, Iloc, 'r-', 'LineWidth', 2);
yrange = [min(Iloc), max(Iday)];
plot(xval(idlock)*ones(1, 2), yrange, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), yrange, '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), yrange, '--', 'Color', grey2, 'LineWidth', 2);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('$I_s$', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);
legend('total', 'local', 'Location', 'best'); legend boxoff;    
if saveTrue
    cd(saveFol);
    saveas(gcf, ['incid' namstr], 'fig');
    cd(thisDir);
end

% Prob(R <= 1) across time
figure;
plot(xval, p1imp, 'b-', 'LineWidth', 2);
hold on;
plot(xval, p1tot, 'r-', 'LineWidth', 2);
plot(xval(idlock)*ones(1, 2), [0 1.05], '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), [0 1.05], '--', 'Color', grey2, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), [0 1.05], '--', 'Color', grey2, 'LineWidth', 2);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('P$(R_s \leq 1)$', 'FontSize', 18); ylim([0 1.05]); 
xlabel('$s$ (days)', 'FontSize', 18);
legend('local', 'total', 'Location', 'best'); legend boxoff;    
if saveTrue
    cd(saveFol);
    saveas(gcf, ['pRless1' namstr], 'fig');
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



