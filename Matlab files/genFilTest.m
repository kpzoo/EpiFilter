% Test EpiFilter vs APE on general scenarios
clearvars; clc; close all; tic;

% Assumptions and notes
% - discrete filter and APE used on same data
% - various epidemic scenarios examined in composite

% Folders for saving data
thisDir = cd; saveFol = 'Results/General';
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Setup simulated epidemics

% Choose scenarios
scenNo = [2 3 5 7]; len = length(scenNo);
%scenNo = [1 4 6 8]; len = length(scenNo);
% SI distribution and if saving data
distNo = 2; saveTrue = 1; 

% Define number of days to simulate
nday0 = 201; tday0 = 1:nday0;
% Epidemic data generated
Iday = cell(1, len); Rtrue = Iday; Lam = Iday;
% No. days (should be same)
nday = zeros(1, len);

for i = 1:len
    % Simulate epidemic scenarios and truncate
    Iwarn = 1; % ensure no warnings
    while Iwarn
        [Iday{i}, Lam{i}, Rtrue{i}, tday, Iwarn, distvals] = epiSimScen(scenNo(i), tday0, nday0, distNo);
    end
    if Iwarn
        warning('Sequences of zero incidence');
    end
    % Truncated observation period
    nday(i) = length(tday);
end

% Check all curves same length
if length(unique(nday)) ~= 1
    error('Incidence curves of inconsistent length');
else
    nday = nday(1);
end
% Saving data and figs name
scenStr = num2str(scenNo); scenStr = scenStr(~isspace(scenStr));
namstr = [num2str(nday) '_' num2str(len) '_' num2str(distNo) '_' scenStr];

%% Recursive filter for each scenario

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);

% Uniform prior over grid of size m
m = 2000; p0 = (1/m)*ones(1, m);
% Estimates for every scenario
Rmed = cell(1, len); Rlow = Rmed; Rhigh = Rmed; Rmean = Rmed;

% Run main filter over every eta
for i = 1:len
    % EpiFilter estimates
    [Rmed{i}, Rlow{i}, Rhigh{i}, Rmean{i}, ~] = runEpiFilter(Rmin, Rmax, m, eta, nday, p0, Lam{i}, Iday{i});
end

%% APE estimate for each scenario

% Estimates for every scenario
Ropt = cell(1, len); Rconf = Ropt; 
% Best windows and ape values
kbest = Ropt; ape = Ropt;

for i = 1:len
    % APEestim estimates
    [~, kbest{i}, Ropt{i}, ~, Rconf{i}, ~, ape{i}] = apeEstim(nday, Iday{i}, Lam{i});
end

%% Plotting and saving

% Colours for plotting
cols = {'r', grey2, 'g', 'b'};
% Plot from second time onwards
tplt = tday(2:end); 

% Four panels of plots
figure;
for j = 1:2:8
   subplot(4, 2, j);
   % APE on left plots
   i = find(j == 1:2:8);
   plot(tplt, Rtrue{i}(2:end), 'k--', 'LineWidth', 2);
   hold on;
   plotCIRaw(tplt', Ropt{i}', Rconf{i}(1, :)', Rconf{i}(2, :)', 'r'); 
   grid off; box off; hold off;
   ylabel(['$\tilde{R}_{\tau(s)} | k = $ ' num2str(kbest{i}(1)+1)], 'FontSize', 18);
   if i == 4
        xlabel('$s$ (days) [APEestim]', 'FontSize', 18);
   end
   xlim([tday(5) tday(end)]); ylim([0 4]);
end
for j = 2:2:8
   subplot(4, 2, j);
   % Filter on right plots
   i = find(j == 2:2:8);
   plot(tplt, Rtrue{i}(2:end), 'k--', 'LineWidth', 2);
   hold on;
   plotCIRaw(tplt', Rmean{i}(2:end)', Rlow{i}(2:end)', Rhigh{i}(2:end)', 'r');
   grid off; box off; hold off;
   ylabel(['$\tilde{R}_s | \eta = $ ' num2str(eta)], 'FontSize', 18);
   if i == 4
        xlabel('$s$ (days) [EpiFilter]', 'FontSize', 18);
   end
   xlim([tday(5) tday(end)]); ylim([0 4]);
end
if saveTrue
    cd(saveFol);
    saveas(gcf, ['compMeth_' namstr], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
if saveTrue
    cd(saveFol);
    save(['genComp' namstr '.mat']);
    cd(thisDir);
end



