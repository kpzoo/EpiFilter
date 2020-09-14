% Write CSVs from empirical studies

% Filtered and smoothed R estimates
csvwrite('RmeanS.csv', Rmean);
csvwrite('RmeanF.csv', RmeanF);
csvwrite('RlowS.csv', Rlow);
csvwrite('RlowF.csv', RlowF);
csvwrite('RhighS.csv', Rhigh);
csvwrite('RhighF.csv', RhighF);

% Filtered and smoothed I predictions
csvwrite('predS.csv', predS);
csvwrite('predF.csv', predF);
csvwrite('predIntS.csv', predIntS);
csvwrite('predIntF.csv', predIntF);