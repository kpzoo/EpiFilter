% Test mechanics of EpiFilter but with varying grid and noise
clearvars; clc; close all; tic;

% Assumptions and notes
% - discrete replicator/Snyder method employed
% - test the grid resolution and state model noise

% Folders for saving data
thisDir = cd; saveFol = 'Results/Mechanics';
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Setup simulated epidemic

% Define number of days to simulate
nday0 = 201; tday0 = 1:nday0;
% Choose scenario, generation time distribution
scenNo = 5; distNo = 2;

% Simulate epidemic scenarios and truncate
Iwarn = 1; % ensure no warnings
while Iwarn
    [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimScen(scenNo, tday0, nday0, distNo);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Truncated observation period
nday = length(tday);

% Saving data and figs
saveTrue = 1; 
namstr = [num2str(nday) '_' num2str(scenNo) '_' num2str(distNo) '_'];

% Initial parameters for EpiFilter
Rmin = 0.01; Rmax = 10; % grid of R

%% Examine various numbers of points 

% Grid point range
m = [20 50 100]; lenm = length(m);
% Noise level for state space
etam = 0.1; disp(['Eta = ' num2str(etam)]);

% Estimates for every m
Rmedm = cell(1, lenm); Rlowm = Rmedm; Rhighm = Rmedm; Rmeanm = Rmedm;

% Run main filter over every m
for i = 1:lenm
    % Uniform prior
    p0 = (1/m(i))*ones(1, m(i));
    % EpiFilter estimates
    [Rmedm{i}, Rlowm{i}, Rhighm{i}, Rmeanm{i}, ~] = runEpiFilter(Rmin, Rmax, m(i), etam, nday, p0, Lam, Iday);
end

%% Examine various noise intensities

% Noise level range
eta = [1/2 1/10 1/20]; lene = length(eta);
% Grid of state space
me = 1000; disp(['m = ' num2str(me)]);
% Uniform prior
p0 = (1/me)*ones(1, me);

% Estimates for every eta
Rmede = cell(1, lene); Rlowe = Rmede; Rhighe = Rmede; Rmeane = Rmede;

% Run main filter over every eta
for i = 1:lene
    % EpiFilter estimates
    [Rmede{i}, Rlowe{i}, Rhighe{i}, Rmeane{i}, ~] = runEpiFilter(Rmin, Rmax, me, eta(i), nday, p0, Lam, Iday);
end

%% Plotting and saving

% Colours for plotting
cols = {'r', grey2, 'b'};

% Estimates with increasing grid resolution
figure;
plot(tday, Rtrue, 'k--', 'LineWidth', 2);
hold on;
for i = 1:lenm
    %plotCIRaw(tday', Rmedm{i}', Rlowm{i}', Rhighm{i}', cols{i});
    stairs(tday', Rmeanm{i}', 'Color', cols{i}, 'LineWidth', 2);
end
grid off; box off; hold off;
ylabel(['$\rilde{R}_{s} | \eta = $ ' num2str(etam)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(5) tday(end)]); ylim([0 5]);

% Estimates with increasing grid resolution
figure;
plot(tday, Rtrue, 'k--', 'LineWidth', 2);
hold on;
for i = 1:lenm
    %plotCIRaw(tday', Rmede{i}', Rlowe{i}', Rhighe{i}', cols{i});
    stairs(tday', Rmeane{i}', 'Color', cols{i}, 'LineWidth', 2);
end
grid off; box off; hold off;
ylabel(['$\tilde{R}_{s} | m = $ ' num2str(me)], 'FontSize', 18)
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(5) tday(end)]); ylim([0 5]);

% Combine both tests into a 6 panel figure
figure;
for j = 1:2:6
    % Grid tests
    subplot(3, 2, j);
    i = find(j == 1:2:6);
    plot(tday, Rtrue, 'k--', 'LineWidth', 2);
    hold on;
    plotCIRaw(tday', Rmeanm{i}', Rlowm{i}', Rhighm{i}', 'r');
    grid off; box off; hold off;
    ylabel(['$\tilde{R}_{s} | m = $ ' num2str(m(i))], 'FontSize', 18)
    if i == 3
        xlabel(['$s$ (days) $| \eta = $ ' num2str(etam)], 'FontSize', 18);
    end
    xlim([tday(5) tday(end)]); ylim([0 5]);
end
for j = 2:2:6
    % Grid tests
    subplot(3, 2, j);
    i = find(j == 2:2:6);
    plot(tday, Rtrue, 'k--', 'LineWidth', 2);
    hold on;
    plotCIRaw(tday', Rmeane{i}', Rlowe{i}', Rhighe{i}', 'r');
    grid off; box off; hold off;
    ylabel(['$\tilde{R}_{s} | \eta = $ ' num2str(eta(i))], 'FontSize', 18)
    if i == 3
        xlabel(['$s$ (days) $| m = $ ' num2str(me)], 'FontSize', 18);
    end
    xlim([tday(5) tday(end)]); ylim([0 5]);
end
if saveTrue
    cd(saveFol);
    saveas(gcf, ['testMech_' namstr], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['mechanics' namstr '.mat']);
    cd(thisDir);
end



