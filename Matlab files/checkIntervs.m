% All interventions and R before and after as well as Z
clc; close all; 

%% Load NZ data and process

load('/Users/kp10/Desktop/Imperial/2020/Term 2/Journals/Elimination paper/figures/proc_New Zealand.mat')

% NZ intervention dates
intervs = {'19-03-20', '26-03-20', '14-05-20', '09-06-20', '12-08-20'};

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

% Extract smoothed estimates of R
Rsmooth = [Rimp.low(:,2) Rimp.med(:,2) Rimp.high(:,2)];
Rnaive = [Rtot.low(:,2) Rtot.med(:,2) Rtot.high(:,2)];

% Extract smoothed estimates of Z
Zsmooth = 100*zsloc'; Znaive = 100*zstot';

% Get times 1-2 weeks after interventions
idwk1 = idint + 7; idwk2 = idint + 14;

% R and Z at intervention times
Rbef = Rsmooth(idint, :); Raf1 =  Rsmooth(idwk1, :); Raf2 =  Rsmooth(idwk2, :);
Zbef = Zsmooth(idint-1); Zaf1 = Zsmooth(idwk1-1); Zaf2 = Zsmooth(idwk2-1);

%% Load HK data and process

load('/Users/kp10/Desktop/Imperial/2020/Term 2/Journals/Elimination paper/figures/proc_Hong Kong.mat')

% HK intervention dates
intervs = {'25-01-20', '25-03-20', '29-03-20', '04-05-20', '05-05-20', '27-05-20', '04-07-20', '13-07-20', '19-07-20', '23-11-20'};

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

% Extract smoothed estimates of R
Rsmooth = [Rimp.low(:,2) Rimp.med(:,2) Rimp.high(:,2)];
Rnaive = [Rtot.low(:,2) Rtot.med(:,2) Rtot.high(:,2)];

% Extract smoothed estimates of Z
Zsmooth = 100*zsloc'; Znaive = 100*zstot';

% Get times 1-2 weeks after interventions
idwk1 = idint + 7; idwk2 = idint + 14;
% Fix if too far in time
idwk1(idwk1 > nday) = nday; idwk2(idwk2 > nday) = nday;

% R and Z at intervention times
Rbef = Rsmooth(idint, :); Raf1 =  Rsmooth(idwk1, :); Raf2 =  Rsmooth(idwk2, :);
Zbef = Zsmooth(idint-1); Zaf1 = Zsmooth(idwk1-1); Zaf2 = Zsmooth(idwk2-1);

%% Load Vic data and process

load('/Users/kp10/Desktop/Imperial/2020/Term 2/Journals/Elimination paper/figures/proc_Vic.mat')

% Victoria intervention dates
intervs = {'16-03-20', '20-03-20', '02-05-20', '14-05-20', '30-06-20', '02-08-20', '18-10-20', '22-11-20', '23-11-20'};

% Reformat intervention dates
intervs = datetime(intervs, 'InputFormat', 'dd-MM-yy');
nInts = length(intervs); idint = zeros(1, nInts);
% Get ids of interventions
for i = 1:nInts
    idint(i) = find(intervs(i) == tdates);
end

% Extract smoothed estimates of R
Rsmooth = [Rimp.low(:,2) Rimp.med(:,2) Rimp.high(:,2)];
Rnaive = [Rtot.low(:,2) Rtot.med(:,2) Rtot.high(:,2)];

% Extract smoothed estimates of Z
Zsmooth = 100*zsloc'; Znaive = 100*zstot';

% Get times 1-2 weeks after interventions
idwk1 = idint + 7; idwk2 = idint + 14;
% Fix if too far in time
idwk1(idwk1 > nday) = nday; idwk2(idwk2 > nday) = nday;

% R and Z at intervention times
Rbef = Rsmooth(idint, :); Raf1 =  Rsmooth(idwk1, :); Raf2 =  Rsmooth(idwk2, :);
Zbef = Zsmooth(idint-1); Zaf1 = Zsmooth(idwk1-1); Zaf2 = Zsmooth(idwk2-1);


