% Plot confidence intervals (asymmetric)
function [l, p] = plotCI2(t, y, e1, e2, colstr, wantStairs)

% Assumptions and notes
% - allows for stairs (piecewise-constant) plots
% - uses boundedline package (assumed on path)
% - all inputs are column vectors, colstr = 'b' for example

% Parse into piecewise-constant form
if wantStairs
    [~, e1] = stairs(t, e1);
    [~, e2] = stairs(t, e2);
    [t, y] = stairs(t, y);
end

% Main CI and median (or mean) plot
[l,p] = boundedline(t, y, [e1 e2], ['-' colstr]);

% Median line width and contrast of CI
l.LineWidth = 2; p.FaceAlpha = 0.5;
l.LineWidth = 2; p.LineStyle = 'none';
% outlinebounds(l,p);