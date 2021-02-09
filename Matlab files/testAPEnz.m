% Examine EpiEstim at various intervention times
clc; close all;

% Assumptions and notes
% - assume processed data from EpiFilter already loaded as .mat
% - distinguishes local vs imported cases

% Directory and if saving
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% NZ intervention dates
intervs = {'19-03-20','26-03-20', '14-05-20', '09-06-20', '12-08-20'};

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

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
    [pred{i}, predInt{i}, prob{i}, R{i}, RInt{i}, wins(i)] = getNegBinAPE2(ks(i), nday, Iloc, Lam);
    
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
klong = 6; disp(['Long window size: ' num2str(klong+1)]);
% Corresponding estimates and predictions
RoptL = R{klong}; IoptL = pred{klong};
RconfL = RInt{klong}; IconfL = predInt{klong};

%% Figures of individual estimates and predictions

% Set k to window length = k + 1 
ks = ks + 1; kbest = kbest + 1; klong = klong+1;
% Plot from second time onwards
tplt = tday(2:end); Iplt = Iday(2:end);

% % One step ahead predictions and metric
% figure;
% subplot(2, 1, 1);
% hold on;
% plot(tplt, ones(size(tplt)), '--', 'color', 'k', 'linewidth', 2);
% plotCIRaw(tplt', Ropt', Rconf(1, :)', Rconf(2, :)', 'r');
% grid off; box off; hold off;
% ylabel(['$\tilde{R}_s | k^* = $ ' num2str(kbest(1))], 'FontSize', 18);
% xlim([tplt(1) tplt(end)]);
% subplot(2, 1, 2);
% hold on;
% scatter(tplt, Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
% plotCIRaw(tplt', Iopt', Iconf(1, :)', Iconf(2, :)', 'r');
% grid off; box off; hold off;
% xlim([tplt(1) tplt(end)]);
% ylabel(['$\tilde{I}_s | k^* = $ ' num2str(kbest(1))], 'FontSize', 18)
% xlabel('$s$ (days)', 'FontSize', 18);
% 
% % Long window predictions and metric
% figure;
% subplot(2, 1, 1);
% hold on;
% plot(tplt, ones(size(tplt)), '--', 'color', 'k', 'linewidth', 2);
% plotCIRaw(tplt', RoptL', RconfL(1, :)', RconfL(2, :)', 'r');
% grid off; box off; hold off;
% ylabel(['$\tilde{R}_s | k^* = $ ' num2str(klong)], 'FontSize', 18);
% xlim([tplt(1) tplt(end)]);
% subplot(2, 1, 2);
% hold on;
% scatter(tplt, Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor', grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
% plotCIRaw(tplt', IoptL', IconfL(1, :)', IconfL(2, :)', 'r');
% grid off; box off; hold off;
% xlim([tplt(1) tplt(end)]);
% ylabel(['$\tilde{I}_s | k^* = $ ' num2str(klong)], 'FontSize', 18)
% xlabel('$s$ (days)', 'FontSize', 18);


%% Publishable figure

% Bayesian smoothed R 
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
ax(1) = subplot(3, 2, 1);
hold on; 
plotCIRaw(xval', Rimp.med(:, 2), Rimp.low(:, 2), Rimp.high(:, 2), 'r');
plot(pltdate, Ryrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (local)', 'FontSize', 18)
ylabel('EpiFilter $R_s$', 'FontSize', 20); ylim([0 3]);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% Predicted causal incidence
ax(2) = subplot(3, 2, 2);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Iimp.med(:, 2), Iimp.low(:, 2), Iimp.high(:, 2), 'b');
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (local)', 'FontSize', 18);
ylabel('EpiFilter $I_s$', 'FontSize', 20);

% APEestim 
ax(3) = subplot(3, 2, 3);
hold on; 
plotCIRaw(xval(2:end)', Ropt', Rconf(1, :)', Rconf(2, :)', 'r');
plot(pltdate, Ryrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (local)', 'FontSize', 18)
ylabel('APEestim $R_s$', 'FontSize', 20); ylim([0 3]);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% Predicted causal incidence
ax(4) = subplot(3, 2, 4);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', Iopt', Iconf(1, :)', Iconf(2, :)', 'b');
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (local)', 'FontSize', 18);
ylabel('APEestim $I_s$', 'FontSize', 20);

% EpiEstim 
ax(5) = subplot(3, 2, 5);
hold on; 
plotCIRaw(xval(2:end)', RoptL', RconfL(1, :)', RconfL(2, :)', 'r');
plot(pltdate, Ryrange, '-', 'Color', grey1, 'LineWidth', 2);
plot(xval, ones(size(tday)), 'k--', 'LineWidth', 1);
grid off; box off; hold off;
%ylabel('$\hat{R}_s$ (local)', 'FontSize', 18)
ylabel('EpiEstim $R_s$', 'FontSize', 20); ylim([0 3]);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
xlabel('time, $s$ (days)', 'FontSize', 18);
% Predicted causal incidence
ax(6) = subplot(3, 2, 6);
scatter(xval(2:end), Iloc(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(xval(2:end)', IoptL', IconfL(1, :)', IconfL(2, :)', 'b');
plot(pltdate, Iyrange, '-', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off;
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
%ylabel('$\hat{I}_s$ (local)', 'FontSize', 18);
ylabel('EpiEstim $I_s$', 'FontSize', 20);
xlabel('time, $s$ (days)', 'FontSize', 18);

if saveTrue
    cd(saveFol);
    saveas(gcf, ['compAll' namstr 'die'], 'fig');
    cd(thisDir);
end

