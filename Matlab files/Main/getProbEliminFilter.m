% Probability an epidemic is eliminated with imports using EpiFilter
function [z0, z1] = getProbEliminFilter(I, Itot, Rgrid, m, eta, distvals, p0)

% Assumptions and notes
% - I can be local or total number of cases (if can't distinguish)
% - Itot is always total number of cases (with imports)
% - append with pseudo-data of 0s to compute
% - only input data up to point want to check elimination
% - posterior for R based on EpiFilter and over Rgrid

% Define a zero look-ahead sequence
Iz = zeros(1, 150); 

% Append epidemic curve with pseudo-data
Icurr = [I Iz]; ncurr = length(Icurr);
% Total cases (local and imported)
Ilamcurr = [Itot Iz];

% Range of index time points
idz = length(I) + 1; ir = idz:ncurr-1; 

% Compute serial interval
serial = serialDistrTypes(ncurr, distvals);
Pomega = serial(1/distvals.omega);

% Compute successive total infectiousness for this I
Lcurr = zeros(1, ncurr);
for i = 2:ncurr
    % Relevant part of SI: Pomega(1:i-1))
    Lcurr(i) = sum(Ilamcurr(i-1:-1:1).*Pomega(1:i-1));
end

% Run EpiFilter approach on relevant incidence and total infectiousness
[~, ~, ~,  RmeanF, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, ncurr, p0, Lcurr, Icurr);
[~, ~, ~,  RmeanS, qR] = runEpiSmoother(Rgrid, m, ncurr, pR, pRup, pstate);

% Estimate probability terms for elimination
zseq0 = exp(-Lcurr(ir+1).*RmeanF(ir));
zseq1 = exp(-Lcurr(ir+1).*RmeanS(ir));

% Estimated probability of elimination at this time
z0 = prod(zseq0); z1 = prod(zseq1);

