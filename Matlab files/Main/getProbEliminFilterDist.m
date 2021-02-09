% Probability an epidemic is eliminated with imports using EpiFilter
function [z0, z1, ztest0, ztest1] = getProbEliminFilterDist(I, Itot, Rgrid, m, eta, distvals, p0)

% Assumptions and notes
% - instead of mean uses the entire distribution
% - I can be local or total number of cases (if can't distinguish)
% - Itot is always total number of cases (with imports)
% - append with pseudo-data of 0s to compute
% - only input data up to point want to check elimination
% - posterior for R based on EpiFilter and over Rgrid

% Define a zero look-ahead sequence
Iz = zeros(1, 100); lz = length(Iz) - 1;

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
ztest0 = prod(exp(-Lcurr(ir+1).*RmeanF(ir)));
ztest1 = prod(exp(-Lcurr(ir+1).*RmeanS(ir)));

% Sequences of probabilities across zeros
zseq0 = zeros(1, lz); zseq1 = zseq0; ix = 1;
for i = ir
    % Lam and elimin prob at this time
    L = Lcurr(i+1); zL = exp(-L*Rgrid);
    % Integrate over grid
    zseq0(ix) = zL*pR(i, :)'; zseq1(ix) = zL*qR(i, :)';
    ix = ix + 1;
end

% Estimated probability of elimination at this time
z0 = prod(zseq0); z1 = prod(zseq1);

