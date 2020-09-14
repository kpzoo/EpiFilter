% Just some default paths
function [grey1, grey2, cmap] = defaultSet(nCol)

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 18);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);
set(0,'DefaultAxesTitleFontWeight','normal');
% New colour choices
cmap = linspecer(nCol, 'sequential');