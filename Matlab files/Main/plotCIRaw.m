% Plot confidence intervals (asymmetric)
function [l, p] = plotCIRaw(t, y, y1, y2, colstr)

% Assumptions and notes
% - uses boundedline package (assumed on path)
% - all inputs are column vectors, colstr = 'b' for example

% Differences from central measure
e1 = y - y1; e2 = y2 - y;

% Main CI and median (or mean) plot
[l,p] = boundedline(t, y, [e1 e2], ['-' colstr]);

% Median line width and contrast of CI
l.LineWidth = 2; p.FaceAlpha = 0.5; p.LineStyle = 'none';
%outlinebounds(l,p);