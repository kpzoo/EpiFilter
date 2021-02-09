% R-Z framework for assessing local transmission and elimination
clearvars; clc; close all; tic;

% Assumptions and notes
% - works with Australian state data
% - includes elimination times and imported vs local
% - provides several serial interval options
% - allows addition of several intervention dates
% - adds publishable figure with incidence, R and Z

% Directory and if saving
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);
grey2 = [0.65 0.65 0.65];

%% Extract empirical data from WHO and setup renewal model

% Possible countries of interest
states = {'Vic', 'SA'};
% Decide country to investigate
scenNo = 1; scenNam = states{scenNo};
disp(['Examining data from state ' scenNam]);

% Identifier for saving
namstr = ['_' scenNam];
% Folders for saving and loading
saveFol = ['Results/' scenNam '/'];
loadFol =  ['Data/' scenNam '/'];

% Data from country of interest
cd(loadFol);

switch scenNo
    case 1
        % Victoria intervention dates
        intervs = {'16-03-20', '30-06-20', '02-08-20', '18-10-20', '23-11-20', '26-11-20'};
        % Start and stop dates
        start = datetime('23-02-20', 'InputFormat', 'dd-MM-yy');
        stop = datetime('27-11-20', 'InputFormat', 'dd-MM-yy');
        
        % Define data-type i.e. histogram or time-series
        dataType = 1;
        switch(dataType)
            case 1
                disp('Using case histograms');
                
                % File with Victoria epidemic curves
                file = dir('VIC.csv'); file = file.name;
                % Data for Victoria state overseas vs local
                data = readtable(file);
                
                % Possible imported and local cases
                Iintro = data.Overseas;
                Iloc = data.KnownLocal + data.UnknownLocal_Community_;
                Iday = Iloc + Iintro;
                
                % Dates of interest
                dates1 = data.Var1; nlen1 = length(dates1);
                dates = cell(1, nlen1);
                for i = 1:nlen1
                    % Add year and reformat
                    dates{i} = [dates1{i} '/20'];
                end
                dates = datetime(dates, 'InputFormat', 'dd/MM/yy');
                
                % Get start and stop and ensure have all in-between
                id1 = find(dates == start); id2 = find(dates == stop);
                if length(id1:id2) ~= length(start:stop)
                    error('Dates not aligned, some might be missing');
                end
                
                % Truncate to range of interest
                Iloc = Iloc(id1:id2)'; Iintro = Iintro(id1:id2)';
                Iday = Iday(id1:id2)'; dates = dates(id1:id2)';
                
            case 2
                disp('Using case time-series');
                
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
                
        end
        
    case 2
        % File with South Australian data
end

% Length of time-series
tdates = start:stop;
nday = length(tdates); tday = 1:nday;
% Numeric value associated with date, 2000 is pivot
xval = datenum(tdates);

cd(thisDir);

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

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

%% Figures and data storage

% All dates as vertical lines
pltdate = kron(xval(idint), ones(1, 2));
pltdate = reshape(pltdate, [2 nInts]);
% Y plotting for Z
Zyrange = 100*[0, 1.05]; Zyrange = reshape(repmat(Zyrange, [1, nInts]), [2 nInts]);

% Y plotting range for R
Ryrange = [min(Rimp.low(:, 2)), max(Rtot.high(:, 2))];
RyrangeF = [min(Rimp.low(:, 1)), max(Rtot.high(:, 1))];
Ryrange = reshape(repmat(Ryrange, [1, nInts]), [2 nInts]);
RyrangeF = reshape(repmat(RyrangeF, [1, nInts]), [2 nInts]);

% Y plotting range for I
Iyrange = [min(Iimp.low(:, 2)), max(Itot.high(:, 2))];
Iyrange = reshape(repmat(Iyrange, [1, nInts]), [2 nInts]);
IyrangeF = [min(Iimp.low(:, 1)), max(Itot.high(:, 1))];
IyrangeF = reshape(repmat(IyrangeF, [1, nInts]), [2 nInts]);


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
plot(pltdate, Zyrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(2:end), 95*ones(size(xval(2:end))), 'k--', 'LineWidth', 1);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('smoothed $Z_s$', 'FontSize', 18); ylim(100*[0 1.05]);
xlabel('time, $s$ (days)', 'FontSize', 18);
legend('naive', 'local', 'Location', 'best'); legend boxoff;
subplot(1, 2, 2);
stairs(xval(2:end), 100*zftot, 'Color', 'g', 'LineWidth', 2);
hold on;
stairs(xval(2:end), 100*zfloc, 'Color', [0.07,0.62,1.00], 'LineWidth', 2);
plot(pltdate, Zyrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval(2:end), 95*ones(size(xval(2:end))), 'k--', 'LineWidth', 1);
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
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('$I_s$', 'FontSize', 18);
xlabel('$s$ (days)', 'FontSize', 18);
legend('total', 'local', 'Location', 'best'); legend boxoff;

% Prob(R <= 1) across time
figure;
plot(xval, p1imp, 'b-', 'LineWidth', 2);
hold on;
plot(xval, p1tot, 'r-', 'LineWidth', 2);
plot(pltdate, Zyrange/100, '-', 'Color', grey1, 'LineWidth', 2);
hold off; box off; grid off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('P$(R_s \leq 1)$', 'FontSize', 18); ylim([0 1.05]);
xlabel('$s$ (days)', 'FontSize', 18);
legend('local', 'total', 'Location', 'best'); legend boxoff;

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    clear('pRup', 'pstate');
    save(['proc' namstr '.mat']);
    cd(thisDir);
end


%% Publishable figure

% Colours and proportion of imports
blu = [0.07,0.62,1.00]; blu2 = [0.06 1 1];
Iprop = Iintro./Iday; Iprop(isnan(Iprop)) = 0;

% Main figure - R(t), z(t) and incidence curve
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);

subplot(3, 1, 1);
% Simple incidence curve
h = area(xval, [Iloc' Iintro']);
h(1).FaceColor = 'r'; h(1).FaceAlpha = 0.6; h(1).EdgeColor = grey1;
h(2).FaceColor = grey1; h(2).FaceAlpha = 0.6; h(2).EdgeColor = grey1;
legend(h, 'local', 'imported', 'Box', 'off', 'autoupdate', 'off');
yrange = [0 max(Iday)+10]; ylim([0 yrange(2)+10]);
hold on;
plot(pltdate, Iyrange, '-', 'Color', grey2, 'LineWidth', 2);hold off; box off; grid off;
ylabel('case incidence', 'FontSize', 18);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
h = gca; h.FontSize = 20;

subplot(3, 1, 2:3);
% Smoothed R estimate with account for imports
yyaxis left; h = gca; h.YColor = 'r';
plot(xval, ones(size(tday)), 'r--', 'LineWidth', 1);
hold on;
plotCIRaw(xval', Rimp.med(:, 2), Rimp.low(:, 2), Rimp.high(:, 2), 'r');
%yrange =  [min(Rimp.low(:, 2)), max(Rimp.high(:, 2))];
yrange = [0, 2.5];
plot(pltdate, Ryrange, '-', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off;
ylabel('local reproduction no. ($R$)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
xlabel('time (days)', 'FontSize', 18);
% Probability of elimination (smoothed) with account for imports
yyaxis right; h = gca; h.YColor = blu;
stairs(xval(2:end), 100*zsloc, 'Color', blu, 'LineWidth', 2);
hold on;
plot(xval(2:end), 95*ones(size(xval(2:end))), '--', 'Color', blu, 'LineWidth', 1);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('\% confidence in elimination ($Z$)', 'FontSize', 18); ylim([0 105]);
grid off; box off; hold off;
h = gca; h.FontSize = 20;

if saveTrue
    cd(saveFol);
    saveas(gcf, ['comb' namstr], 'fig');
    cd(thisDir);
end
