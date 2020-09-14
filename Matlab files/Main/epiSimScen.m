% Simulate epidemic via renewal model
function [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimScen(scenNo, tday, nday, distNo)

% Assumptions and notes
% - removes 'burn-in' of first 20 days, epidemic size < 100
% - various R trajectories adns SI distributions specified

% Possible scenarios available - must match calling function
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};
disp(['True R scenario: ' scenNam{scenNo}]);

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

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo}; disp(['True SI scenario: ' distChoice]);
% Set serial interval parameters
distvals.type = distNo; 

% Hyerparameters of serial distribution
switch(distNo)
    case 1
        % Geometric distribution - no parameter
        distvals.pm = []; distvals.omega = 15.3;
    case 2
        % Gamma distribution - shape parameter
        %distvals.pm = (1/0.65)^2; distvals.omega = 6.5;
        distvals.pm = 2.7066; distvals.omega = 15.3;
    case 3
        % Delta distribution - odd window around mean
        distvals.pm = 7; distvals.omega = 15.3;
    case 4
        % Two Gamma distributions for flare-up
        distvals.pm = 45; distvals.omega = 15.3;
end

% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Daily incidence and infectiousness
Iday = zeros(size(nday)); Lam = Iday; 
% Initialise epidemic and warning
Iday(1) = 10; Iwarn = 0; 

% Iteratively generate renewal epidemic
for i = 2:nday
    % Relevant part of serial distribution
    Pomegat = Pomega(1:i-1);
    % Total infectiousness
    Lam(i) = sum(Iday(i-1:-1:1).*Pomegat);    
    % Renewal incidence
    Iday(i) = poissrnd(Rtrue(i)*Lam(i));
end

% Remove start-up 20 days
idz = 20:nday; tday = idz;
% Adjusted vectors - including tday
Iday = Iday(idz); Rtrue = Rtrue(idz); Lam = Lam(idz);

% Remove small epidemics
if sum(Iday) < 100
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end
