% Publishable figures assuming output of localCOVIDElim.m
close all; clc;

% HK intervention dates
intervs = {'25-01-20', '25-03-20', '29-03-20', '05-05-20', '27-05-20', '04-07-20', '13-07-20', '20-07-20', '23-11-20'};

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
plot(pltdate, Iyrange, '-', 'Color', grey2, 'LineWidth', 2);
hold off; box off; grid off;
ylabel('case incidence', 'FontSize', 18);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
h = gca; h.FontSize = 22;

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
plot(xval(2:end), 50*ones(size(xval(2:end))), '--', 'Color', blu, 'LineWidth', 1);
xlim([xval(2) xval(end)]); datetick('x','dd-mm', 'keeplimits');
ylabel('\% confidence in elimination ($Z$)', 'FontSize', 18); ylim([0 105]);
grid off; box off; hold off;
h = gca; h.FontSize = 22;

if saveTrue
    cd(saveFol);
    saveas(gcf, ['comb' namstr], 'fig');
    cd(thisDir);
end
