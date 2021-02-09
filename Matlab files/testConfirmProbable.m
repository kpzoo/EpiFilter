% Quick comparison code
clearvars; clc; close all; 
thisDir = cd; saveTrue = 1;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

% Load confirmed cases
load('/Users/kp10/Desktop/Imperial/2020/Term 1/Code/EpiFilter NZ/Matlab files/Results/COVID/new neil serial/proc_New Zealand.mat');

% Incidence
Iday1 = Iday; Iintro1 = Iintro; Iloc1 = Iloc;

% R numbers (medians)
Rloc1 = Rimp.med(:, 2); Rtot1 =  Rtot.med(:, 2);

% Z numbers
zsloc1 = zsloc; zfloc1 = zfloc; 
zstot1 = zstot; zftot1 = zftot;

% Load confirmed and probable cases
load('/Users/kp10/Desktop/Imperial/2020/Term 1/Code/EpiFilter NZ/Matlab files/Results/confirm_probable/proc_New Zealand.mat');

% Incidence
Iday2 = Iday; Iintro2 = Iintro; Iloc2 = Iloc;

% R numbers (medians)
Rloc2 = Rimp.med(:, 2); Rtot2 =  Rtot.med(:, 2);

% Z numbers
zsloc2 = zsloc; zfloc2 = zfloc; 
zstot2 = zstot; zftot2 = zftot;


% Declatation times
tdecLoc = [find(zsloc1 >= 0.95, 1, 'first'), find(zsloc2 >= 0.95, 1, 'first')];
tdecTot = [find(zstot1 >= 0.95, 1, 'first'), find(zstot2 >= 0.95, 1, 'first')];

% Compare Z
figure;
subplot(2, 1, 1);
plot(tday(2:end), zsloc1, tday(2:end), zsloc2, 'LineWidth', 2);
grid off; box off;
ylabel('$Z_s$ (local)', 'FontSize', 18); 
subplot(2, 1, 2);
plot(tday(2:end), zstot1, tday(2:end), zstot2, 'LineWidth', 2);
grid off; box off;
ylabel('$Z_s$ (naive)', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);

% Compare I
figure;
subplot(3, 1, 1);
plot(tday, Iloc1, tday, Iloc2, 'LineWidth', 2);
grid off; box off;
ylabel('$L_s$', 'FontSize', 18); 
subplot(3, 1, 2);
plot(tday, Iintro1, tday, Iintro2, 'LineWidth', 2);
grid off; box off;
ylabel('$M_s$', 'FontSize', 18); 
subplot(3, 1, 3);
plot(tday, Iday1, tday, Iday2, 'LineWidth', 2);
grid off; box off;
ylabel('$I_s$', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);

% Compare difference in I
figure;
subplot(3, 1, 1);
plot(tday, Iloc2 - Iloc1, 'LineWidth', 2);
grid off; box off;
ylabel('$\Delta L_s$', 'FontSize', 18); 
subplot(3, 1, 2);
plot(tday, Iintro2 - Iintro2, 'LineWidth', 2);
grid off; box off;
ylabel('$\Delta M_s$', 'FontSize', 18); 
subplot(3, 1, 3);
plot(tday, Iday2 - Iday1, 'LineWidth', 2);
grid off; box off;
ylabel('$\Delta I_s$', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);


% Compare R
figure;
subplot(2, 1, 1);
plot(tday, Rloc1, tday, Rloc2, 'LineWidth', 2);
grid off; box off;
ylabel('$R_s$ (local)', 'FontSize', 18); 
subplot(2, 1, 2);
plot(tday, Rtot1, tday, Rtot2, 'LineWidth', 2);
grid off; box off;
ylabel('$R_s$ (naive)', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);

% Compare Z, R and I 
figure;
subplot(3, 1, 1);
plot(tday, Iday1, tday, Iday2, 'LineWidth', 2);
grid off; box off;
ylabel('$I_s$ (total)', 'FontSize', 18); 
xlim([tday(2) tday(end)]);
subplot(3, 1, 2);
plot(tday, Rloc1, tday, Rloc2, 'LineWidth', 2);
grid off; box off;
ylabel('$R_s$ (local)', 'FontSize', 18); 
xlim([tday(2) tday(end)]);
subplot(3, 1, 3);
plot(tday(2:end), zsloc1, tday(2:end), zsloc2, 'LineWidth', 2);
grid off; box off;
ylabel('$Z_s$ (local)', 'FontSize', 18); 
xlabel('$s$ (days)', 'FontSize', 18);
xlim([tday(2) tday(end)]);

