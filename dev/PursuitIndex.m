function [PIS] = PursuitIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, TI)

DataPath = GetDataPath();

filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.DPI.mat');
Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
Expt = Neuron.Expt;

StartTime = 500;
FinishTime = median([Expt.Trials(:).End] - [Expt.Trials(:).Start])+StartTime;

dxValues = unique([Expt.Trials.dx]);

SpikeCounts = zeros(length([Expt.Trials]),1);
for tr = 1: length([Expt.Trials]), 
    SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
end

if isfield(Expt.Trials,'dfx')
    deltafxy = [Expt.Trials(:).dfx] - [Expt.Trials(:).fx];
else
    deltafxy = [Expt.Trials(:).dfy] - [Expt.Trials(:).fy];
end

if isfield(Expt.Trials,'dfx')
    dfx = [Expt.Trials(:).dfx] - [Expt.Trials(:).fx];
else
    dfx = zeros(1,length(Expt.Trials));
end
if isfield(Expt.Trials,'dfy')
    dfy = [Expt.Trials(:).dfy] - [Expt.Trials(:).fy];
else
    dfy = zeros(1,length(Expt.Trials));
end

deltatheta = atan2(dfy, dfx) + pi/2 - Expt.Stimvals.or .* pi / 180;
%TEST 
disp([' = = = =>> ', Expt.Header.Name, ' - or =' , num2str(Expt.Stimvals.or) ,' TEST = ', num2str(sum((cos(deltatheta)>0) ~= (mod(deltatheta, 2 * pi)==0)))]);
%deltatheta = mod(deltatheta, 2 * pi);
deltatheta = (cos(deltatheta)>0);
    

deltafxy = fix(deltafxy);
Speeds = unique(abs(deltafxy(deltafxy>0)));
disp(Speeds);

C__1 = (deltatheta == 0);
C__2 = (deltatheta ~= 0);

C1_1 = (deltatheta == 0) & (abs(deltafxy)== max(Speeds));
C1_2 = (deltatheta ~= 0) & (abs(deltafxy) ==max(Speeds));

C2_1 = (deltatheta == 0) & (abs(deltafxy)== min(Speeds));
C2_2 = (deltatheta ~= 0) & (abs(deltafxy)== min(Speeds));

C3_1 = (deltatheta == 0) & (abs(deltafxy)==Speeds(round(end/2)));
C3_2 = (deltatheta ~= 0) & (abs(deltafxy)==Speeds(round(end/2)));


C101 = (([Expt.Trials(:).dx] == 0) & (deltatheta == 0) & (abs(deltafxy)== max(Speeds)) );
C102 = (([Expt.Trials(:).dx] == 0) & (deltatheta ~= 0) & (abs(deltafxy)== max(Speeds)) );

C201 = (([Expt.Trials(:).dx] == 0) & (deltatheta == 0) & (abs(deltafxy)== min(Speeds)) );
C202 = (([Expt.Trials(:).dx] == 0) & (deltatheta ~= 0) & (abs(deltafxy)== min(Speeds)));

C301 = (([Expt.Trials(:).dx] == 0) & (deltatheta == 0) & (abs(deltafxy)==Speeds(round(end/2))) );
C302 = (([Expt.Trials(:).dx] == 0) & (deltatheta ~= 0) & (abs(deltafxy)==Speeds(round(end/2))));


C111 = (([Expt.Trials(:).dx] > 0) & (deltatheta == 0) & (abs(deltafxy)== max(Speeds)) );
C112 = (([Expt.Trials(:).dx] > 0) & (deltatheta ~= 0) & (abs(deltafxy)== max(Speeds)) );

C211 = (([Expt.Trials(:).dx] > 0) & (deltatheta == 0) & (abs(deltafxy)== min(Speeds)) );
C212 = (([Expt.Trials(:).dx] > 0) & (deltatheta ~= 0) & (abs(deltafxy)== min(Speeds)));

C311 = (([Expt.Trials(:).dx] > 0) & (deltatheta == 0) & (abs(deltafxy)==Speeds(round(end/2))) );
C312 = (([Expt.Trials(:).dx] > 0) & (deltatheta ~= 0) & (abs(deltafxy)==Speeds(round(end/2))));


C121 = (([Expt.Trials(:).dx] < 0) & (deltatheta == 0) & (abs(deltafxy)== max(Speeds)) );
C122 = (([Expt.Trials(:).dx] < 0) & (deltatheta ~= 0) & (abs(deltafxy)== max(Speeds)) );

C221 = (([Expt.Trials(:).dx] < 0) & (deltatheta == 0) & (abs(deltafxy)== min(Speeds)) );
C222 = (([Expt.Trials(:).dx] < 0) & (deltatheta ~= 0) & (abs(deltafxy)== min(Speeds)));

C321 = (([Expt.Trials(:).dx] < 0) & (deltatheta == 0) & (abs(deltafxy)==Speeds(round(end/2))) );
C322 = (([Expt.Trials(:).dx] < 0) & (deltatheta ~= 0) & (abs(deltafxy)==Speeds(round(end/2))));


PIFast0 = (mean(SpikeCounts(C101)) - mean(SpikeCounts(C102))) ./ (mean(SpikeCounts(C101)) + mean(SpikeCounts(C102)));
PIFast1 = (mean(SpikeCounts(C111)) - mean(SpikeCounts(C112))) ./ (mean(SpikeCounts(C111)) + mean(SpikeCounts(C112)));
PIFast2 = (mean(SpikeCounts(C121)) - mean(SpikeCounts(C122))) ./ (mean(SpikeCounts(C121)) + mean(SpikeCounts(C122)));

PISlow0 = (mean(SpikeCounts(C201)) - mean(SpikeCounts(C202))) ./ (mean(SpikeCounts(C201)) + mean(SpikeCounts(C202)));
PISlow1 = (mean(SpikeCounts(C211)) - mean(SpikeCounts(C212))) ./ (mean(SpikeCounts(C211)) + mean(SpikeCounts(C212)));
PISlow2 = (mean(SpikeCounts(C221)) - mean(SpikeCounts(C222))) ./ (mean(SpikeCounts(C221)) + mean(SpikeCounts(C222)));

PIMid0 = (mean(SpikeCounts(C301)) - mean(SpikeCounts(C302))) ./ (mean(SpikeCounts(C301)) + mean(SpikeCounts(C302)));
PIMid1 = (mean(SpikeCounts(C311)) - mean(SpikeCounts(C312))) ./ (mean(SpikeCounts(C311)) + mean(SpikeCounts(C312)));
PIMid2 = (mean(SpikeCounts(C321)) - mean(SpikeCounts(C322))) ./ (mean(SpikeCounts(C321)) + mean(SpikeCounts(C322)));

PI     = (mean(SpikeCounts(C__1)) - mean(SpikeCounts(C__2))) ./ (mean(SpikeCounts(C__1)) + mean(SpikeCounts(C__2)));
PIFast = (mean(SpikeCounts(C1_1)) - mean(SpikeCounts(C1_2))) ./ (mean(SpikeCounts(C1_1)) + mean(SpikeCounts(C1_2)));
PISlow = (mean(SpikeCounts(C2_1)) - mean(SpikeCounts(C2_2))) ./ (mean(SpikeCounts(C2_1)) + mean(SpikeCounts(C2_2)));
PIMid  = (mean(SpikeCounts(C3_1)) - mean(SpikeCounts(C3_2))) ./ (mean(SpikeCounts(C3_1)) + mean(SpikeCounts(C3_2)));



D1_1 = (([Expt.Trials(:).dx] > 0) & (deltatheta == 0));
D1_2 = (([Expt.Trials(:).dx] > 0) & (deltatheta ~= 0));

D2_1 = (([Expt.Trials(:).dx] < 0) & (deltatheta == 0));
D2_2 = (([Expt.Trials(:).dx] < 0) & (deltatheta ~= 0));

D0_1 = (([Expt.Trials(:).dx] == 0) & (deltatheta == 0));
D0_2 = (([Expt.Trials(:).dx] == 0) & (deltatheta ~= 0));

dx_1 = (mean(SpikeCounts(D1_1)) - mean(SpikeCounts(D1_2))) ./ (mean(SpikeCounts(D1_1)) + mean(SpikeCounts(D1_2)));
dx_2 = (mean(SpikeCounts(D2_1)) - mean(SpikeCounts(D2_2))) ./ (mean(SpikeCounts(D2_1)) + mean(SpikeCounts(D2_2)));
dx_0 = (mean(SpikeCounts(D0_1)) - mean(SpikeCounts(D0_2))) ./ (mean(SpikeCounts(D0_1)) + mean(SpikeCounts(D0_2)));

if(TI<0)
  dx_Pref = -dx_2;
  dx_Null = -dx_1;
else
  dx_Pref = -dx_1;
  dx_Null = -dx_2;
end

dxEffect  = ((mean(SpikeCounts(D1_1)) - mean(SpikeCounts(D2_1))) + (mean(SpikeCounts(D1_2)) - mean(SpikeCounts(D2_2)))) ./ ...
            ((mean(SpikeCounts(D1_1)) + mean(SpikeCounts(D2_1))) + (mean(SpikeCounts(D1_2)) + mean(SpikeCounts(D2_2))));

x = abs([(mean(SpikeCounts(C101)) - mean(SpikeCounts(C102))) , ...
                  (mean(SpikeCounts(C111)) - mean(SpikeCounts(C112))) , ...
                  (mean(SpikeCounts(C121)) - mean(SpikeCounts(C122)))]);
x = x(~isnan(x));
if (isempty(x))
    dxMaxDelta = -1;
else
    dxMaxDelta = max(x) ./ ((FinishTime - StartTime) ./ 10000);
end
%[B, Brate] = CalcPursuitIndex(Expt);

PIS= [-PIFast2, -PIFast0, -PIFast1, ...
      -PISlow2, -PISlow0, -PISlow1, ...
      -PIMid2, -PIMid0, -PIMid1, ...
      PI, ... %10
      -PIFast, -PISlow, -PIMid, ...
      dxEffect, dxMaxDelta, ... % 15
      dx_1, dx_2, dx_0, dx_Pref, dx_Null]; % 20
 
%disp(B);
%disp(PIS);

end










