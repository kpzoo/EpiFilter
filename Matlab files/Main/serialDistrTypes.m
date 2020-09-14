function pdistr = serialDistrTypes(tmax, distvals)

% Assumptions and notes
% - Chooses between discrete distributions on days
% - Insert max days and calculate serial probabilities
% - p must be a parameter, tday an array of integers
% - p is 1/mean of each distribution


switch(distvals.type)
    case 1
        % Geometric distribution over tmax for a given p
        pdistr = @(p) geomDistr(p, 1:tmax);
    case 2
        % Gamma distribution with integer shape (Erlang)
        pdistr = @(p) gammaDistr(p, 1:tmax, distvals.pm);
    case 3 
        % Delta distribution around p
        pdistr = @(p) deltaDistr(p, tmax, distvals.pm);
    case 4
        % Flare-up from secondary transmission (Lee 2019)
        pdistr = @(p) flareDistr(p, 1:tmax, distvals.pm);
        
    otherwise
        disp('No valid distribution specified');
        return;
end


%% Gamma distribution, p is 1/mean, shape param preset
function pr = gammaDistr(p, x, shapePm)

% shapePm is shape parameter of gamma
% default of 20

% Scale parameter based on mean = shapePm*scalePm
scalePm = 1/(p*shapePm); ratePm = 1/scalePm;

% Gamma (Erlang) probabilities
pr = -log(gamma(shapePm)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(x) - ratePm*x;
pr = exp(pr);
%pr = gampdf(k, shapePm, scalePm);

%% Geometric distribution, p is prob success
function pr = geomDistr(p, x)

% Set p to same dimension as x domain
p1 = p;
p = p*ones(size(x));

% Geometric distribution starting at 1 (vs 0)
pr = realpow((1-p), (x-1));
pr = p1*pr;

%% Delta distribution, 1/p is mean
function pr = deltaDistr(p, xmax, nWin)

% nWin is window of mass (must be odd) around mean 
% default of 11

% Distribution defined over all time
pr = zeros(1, xmax);

% No. time points either side of omega
a = round(1/p) - (nWin - 1)/2;
b = round(1/p) + (nWin - 1)/2;

% Probability over this window
pr(a:b) = 1/nWin;

%% Flare-up distribution mixture
function pr = flareDistr(p, x, shapePm)

% shapePm is shape parameter of gamma (default 20)

% Scale parameter based on mean = shapePm*scalePm
scalePm = 1/(p*shapePm); ratePm = 1/scalePm;
% Gamma (Erlang) probabilities
pr = -log(factorial(shapePm-1)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(x) - ratePm*x;
pr = exp(pr);

% Flare-up second distribution - larger mean
flareMean = 40; p2 = 1/flareMean;
scalePm = 1/(p2*shapePm); ratePm = 1/scalePm;
pr2 = -log(factorial(shapePm-1)) + shapePm*log(ratePm) +...
    (shapePm-1)*log(x) - ratePm*x;
pr2 = exp(pr2);

% Combined distribution
pr = pr + pr2;
pr = pr/sum(pr);

