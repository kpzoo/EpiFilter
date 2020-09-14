% Simulate dying epidemic via renewal model
function [Iday, Lam, Rtrue, tday, Iwarn, distvals] = epiSimScenDie(scenNo, nday, distNo, simVals)

% Assumptions and notes
% - removes 'burn-in' of first 20 days, epidemic size < 100
% - various R trajectories adns SI distributions specified

% Possible scenarios available - must match calling function
scenNam = {'control', 'square-wave', 'cascade', 'boom-bust', 'filtered', 'waves', 'noise valley', 'boom-bust-boom'};
disp(['True R scenario: ' scenNam{scenNo}]);

% Parameters for changepoints in R trajectory
Rch = simVals.Rch; tch = simVals.tch;
% Variable for true R
Rtrue = zeros(1, nday);

% Functions for scenarios: R on a daily basis
switch(scenNo)
     case 1
        % Rapidly controlled epidemic
        Rtrue(1:tch) = Rch(1);
        Rtrue(tch+1:end) = Rch(2);
    case 2
        % Rapid control that recovers
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
    case 3
        % Three stage control with fluctuations
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:tch(3)) = Rch(3);
        Rtrue(tch(3)+1:end) = Rch(4);
        Rtrue = Rtrue + 0.3*cosd(2*(1:nday));
        Rtrue(Rtrue <= 0) = 0.1;
    case 4
        % Exponential rise and fall
        trise = 1:tch; tfall = tch+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tch)); Rmax = Rtrue(tch);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.008*(tfall - tch));
    case 5
        % Two stage control with filtered noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
        % Check for noise
        if any(Rtrue < 0)
            Rtrue(Rtrue < 0) = 0;
        end
    case 6
        % Second (sine) wave dynamics
        Rtrue = Rch(1) + Rch(2)*sind(tch(1)*(1:nday));
    case 7
        % Long period of low R between transmission and noise
        Rtrue(1:tch(1)) = Rch(1);
        Rtrue(tch(1)+1:tch(2)) = Rch(2);
        Rtrue(tch(2)+1:end) = Rch(3);
        % Add noise and filter
        Rtrue = Rtrue + normrnd(0.5, 2, [1 nday]);
        m = 20; B = ones(m, 1)/m;
        Rtrue = filter(B, 1, Rtrue);
        % Check for noise
        if any(Rtrue < 0)
            Rtrue(Rtrue < 0) = 0;
        end
    case 8
        % Exponential rise and fall then rise
        tch1 = tch(1); tch2 = tch(2);
        trise = 1:tch1; tfall = tch1+1:tch2; 
        triseSec = tch2+1:nday; % second wave
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.03*(1:tch1)); Rmax = Rtrue(tch1);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.015*(tfall - tch1));
        % Second wave
        Rtrue(triseSec) = Rtrue(tch2)*exp(0.02*(triseSec - tch2));
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
        %distvals.pm = (1/0.65)^2; distvals.omega = 6.5
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
Iday = zeros(1, nday); Lam = Iday; 
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
