% Simulate epidemic via renewal model
function [Iday, Lam, Rtrue, tday, Iwarn] = epiSimAPE2(scenNo, tday, nday, remGap)

% Assumptions and notes
% - warning statements suppressed
% - option to remove a startup sequence of zeros
% - warns if sequence of consecutive zero incidences
% - various R trajectories specified


% Possible scenarios available - must match calling function
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};
disp(['True R scenario: ' scenNam{scenNo}]);
% Warning about incidence zeros
Iwarn = 0; % will be set conditionally

% Functions for scenarios: R on a daily basis
switch(scenNo)
    case 1
        % Constant R
        Rtrue = 1.5*ones(1, nday);
    case 2
        % Sinusoidal R across time
        Rtrue = 1.2 + 0.8*sind(3*tday);
        %Rtrue = 1.5 + 0.5*sind(3*tday);
    case 3
        % Logistic fall in R with time
        R0 = 1.1; Rd = 0.7;
        t1 = 1; t2 = floor(nday/2 + 20);
        Rtrue = R0 + Rd*((1 + exp(-t1*t2))./(1 + exp(-t1*(t2 - tday))));
    case 4
        % Switch point halfway
        tch = ceil(nday/2);
        Rtrue = zeros(1, nday);
        Rtrue(1:tch) = 2.5;
        Rtrue(tch+1:end) = 0.5;
    case 5
        % Exponential rise and fall
        Rtrue = zeros(size(tday)); tchange = floor(nday/4);
        trise = 1:tchange; tfall = tchange+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tchange)); Rmax = Rtrue(tchange);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.008*(tfall - tchange));
    case 6
        % Bottleneck
        Rtrue = zeros(size(tday)); tchange = floor(nday/4);
        Rtrue(1:tchange) = 2;
        Rtrue(tchange+1:2*tchange) = 1.1;
        Rtrue(2*tchange+1:end) = 0.5;
    case 7
        % Two major swings
        Rtrue = zeros(size(tday)); tchange = floor(nday/3);
        Rtrue(1:tchange) = 1;
        Rtrue(tchange+1:2*tchange) = 2;
        Rtrue(2*tchange+1:end) = 0.5;
    case 8
        % White noise with reflection
        Rtrue = abs(normrnd(0.5, 2, [1 nday]));
        % Smooth with m point averager
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
end

% Serial distribution over all tday
distType = 2; %(1 = geom, 2 = erlang)
serial = serialDistrs(nday, distType);
% Single omega controls distribution
omega = 14.2; 
Pomega = serial(1/omega);

% Daily incidence
Iday = zeros(size(nday)); 
% Infectiousness, Poisson rate 
Lam = Iday; rate = Iday;
% Initialise epidemic
Iday(1) = 10;

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = Iday(i-1:-1:1)*Pomegat';
    % Rate for ith day incidence
    rate(i) = Lam(i)*Rtrue(i);
    % Renewal incidence
    Iday(i) = poissrnd(rate(i));
end

% Gaps between non-zero indicence values
zerogaps = diff(find(Iday ~= 0));
% Remove startup sequence of zeros if big
z1 = zerogaps(1);
if z1 > 5 && remGap
    % Update incidence and related vectors
    idz = z1+1:nday;
    % Flag zero incidence regions after startup
    if max(zerogaps(2:end)) > 5
        %warning('Zero incidences beyond startup');
        Iwarn = 1;
    end
else
    % Flag any zero incidence region
    if max(zerogaps) > 5
        %warning('Sequences of zero incidence');
        Iwarn = 1;
    end
    % Un-truncated day set
    idz = 2:nday;
end

% Adjusted vectors - including tday
clear('tday'); tday = idz;
Iday = Iday(idz); 
Rtrue = Rtrue(idz); Lam = Lam(idz);
