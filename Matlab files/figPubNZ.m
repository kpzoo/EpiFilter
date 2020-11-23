% Publishable figures assuming output of nzCOVIDElim.m
close all; clc;

% Colours and proportion of imports
blu = [0.07,0.62,1.00]; blu2 = [0.06 1 1];
Iprop = Iintro./Iday; Iprop(isnan(Iprop)) = 0;

% Main figure - R(t), z(t) and incidence curve
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
% subplot(3, 1, 1);
% % Simple incidence curve with difference
% yyaxis left; h = gca; h.YColor = 'r';
% stairs(xval, Iday, 'Color', blu, 'LineWidth', 2);
% hold on;
% stairs(xval, Iloc, 'r-', 'LineWidth', 2);
% hold off; box off; grid off;
% xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
% ylabel('local cases', 'FontSize', 18);
% ylim([0 max(Iday) + 20]);
% % Proportion of imports
% yyaxis right; h = gca; h.YColor = blu;
% stairs(xval, Iprop, 'Color', blu, 'LineWidth', 1);
% hold on; yrange = [0, 1.1];
% plot(xval(idlock)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
% plot(xval(idrelax)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
% plot(xval(idsecwv)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
% hold off; box off; grid off;
% ylim([0 1.2]); ylabel('imported fraction', 'FontSize', 18);

subplot(3, 1, 1);
% Simple incidence curve 
h = area(xval, [Iloc' Iintro']); 
legend(h, 'local', 'imported', 'Box', 'off', 'autoupdate', 'off');
yrange = [0 max(Iday)+10]; ylim([0 yrange(2)+10]);
hold on;
plot(xval(idlock)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
hold off; box off; grid off;
ylabel('no. cases', 'FontSize', 18);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');

subplot(3, 1, 2:3);
% Smoothed R estimate with account for imports
yyaxis left; h = gca; h.YColor = 'r';
plot(xval, ones(size(tday)), 'r--', 'LineWidth', 1);
hold on; 
plotCIRaw(xval', Rimp.med(:, 2), Rimp.low(:, 2), Rimp.high(:, 2), 'r');
%yrange =  [min(Rimp.low(:, 2)), max(Rimp.high(:, 2))]; 
yrange = [0, 2.5];
plot(xval(idlock)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
plot(xval(idrelax)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
plot(xval(idsecwv)*ones(1, 2), yrange, '-', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off;
ylabel('local reproduction no. ($R$)', 'FontSize', 18)
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
xlabel('time (days)', 'FontSize', 18);
% Probability of elimination (smoothed) with account for imports
yyaxis right; h = gca; h.YColor = blu;
stairs(xval(2:end), zsloc, 'Color', blu, 'LineWidth', 2);
hold on;
plot(xval(2:end), 0.95*ones(size(xval(2:end))), '--', 'Color', blu, 'LineWidth', 1);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('prob. elimination ($Z$)', 'FontSize', 18); ylim([0 1.05]); 
grid off; box off; hold off;