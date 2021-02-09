% Examine early-warning R and Z signals at various intervention times
clc; close all;

% Assumptions and notes
% - assume processed data from EpiFilter already loaded as .mat

% Directory and if saving
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% HK intervention dates
intervs = {'25-01-20', '29-03-20', '05-05-20', '27-05-20', '04-07-20', '13-07-20', '20-07-20', '23-11-20'};
intervs = intervs([3 5 7 8]);

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

% Define indices and portions of epidemic curve
idends = unique(idint);
nends = length(idends); idset = cell(1, nends); 
% Store relevant variables
Iset = idset; Lset = idset; tset = idset; 
Rmedset = idset; Rlowset = idset; Rhighset = idset;
predset = idset; predIntset = idset; xset = idset;

% For each portion compute estimates and predictions
for i = 1:nends
    % Portion considered 
    idset{i} = 1:idends(i);
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
    ylabel('$\hat{I}_s$', 'FontSize', 20);
    xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
    if i == nends
        xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 20);
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
ylabel('$\hat{R}_s$', 'FontSize', 20);
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
xlabel(['$s$ (days) $| \eta = $ ' num2str(eta)], 'FontSize', 20);