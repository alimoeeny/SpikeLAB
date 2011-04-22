% eye calibration

clear, clc
DataPath = GetDataPath();

% % % % Eye Calib
load('../AllEyeCallExpts.mat');
AllNeurons = AllEyeCallExpts;
clear AllEyeCallExpts;
FileType = 'EyeCal';
StimulusType = '';

% % DID
% experiments with valid eyetracking for  both eyes
%BinocularExperiments = {'dae014';'dae015';'dae018';'dae054';'dae149';'dae151';'dae447';'dae4472';'dae449';'dae4492';};



%par
for iN= 1:length(AllNeurons), 
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , num2str(rem(now,1))));
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'); 
    filename(findstr(filename, '..')) = '';
    %disp(filename);
    
    if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', filename), 'file')==2)
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;
        
         lv = zeros(length(Expt.Trials), 2000); lh = zeros(length(Expt.Trials), 2000);
        rv = zeros(length(Expt.Trials), 2000); rh = zeros(length(Expt.Trials), 2000);
        for i = 1: length(Expt.Trials)
            validpoints = Expt.Header.emtimes>0 & Expt.Header.emtimes < 25000;
            lh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,1);
            lv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,3);
            rh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,2);
            rv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,4);
        end
        
    
        scaling(iN) = Expt.Header.emtimes(end) / length(Expt.Header.emtimes) / 10; %Expt.Header.emtimes(end);
        disp([ ' -----------------------------' , num2str(scaling(iN))]);
                        
        R = [Expt.Trials(:).fx] > 0;
        T = [Expt.Trials(:).fy] > 0;
        L = [Expt.Trials(:).fx] < 0;
        B = [Expt.Trials(:).fy] < 0;
        Z = ([Expt.Trials(:).fx] == 0) & ([Expt.Trials(:).fy] == 0);
        
        disp((sum(R) + sum(L) + sum(T) + sum(B) + sum(Z)) / length(Expt.Trials));
        
        Zhr = mean(rh(Z));
        Zvr = mean(rv(Z));
        Zhl = mean(lh(Z));
        Zvl = mean(lv(Z));

        
        Rhr = mean(rh(R)) - Zhr;
        Rvr = mean(rv(R)) - Zvr;
        Rhl = mean(lh(R)) - Zhl;
        Rvl = mean(lv(R)) - Zvl;

        Lhr = mean(rh(L)) - Zhr;
        Lvr = mean(rv(L)) - Zvr;
        Lhl = mean(lh(L)) - Zhl;
        Lvl = mean(lv(L)) - Zvl;
        
        Thr = mean(rh(T)) - Zhr;
        Tvr = mean(rv(T)) - Zvr;
        Thl = mean(lh(T)) - Zhl;
        Tvl = mean(lv(T)) - Zvl;

        Bhr = mean(rh(B)) - Zhr;
        Bvr = mean(rv(B)) - Zvr;
        Bhl = mean(lh(B)) - Zhl;
        Bvl = mean(lv(B)) - Zvl;        
        
        
        figure(1112), clf, hold on,
        %scatter([Rhr, Lhr, Thr, Bhr, Zhr], [Rvr, Lvr, Tvr, Bvr, Zvr], 'r', 'filled');
        %scatter([Rhl, Lhl, Thl, Bhl, Zhl], [Rvl, Lvl, Tvl, Bvl, Zvl], 'b', 'filled');
        scatter([Rhr], [Rvr], 'r', 'filled', 'o');
        scatter([Rhl], [Rvl], 'r', 'filled', 'v');
        scatter([Lhr], [Lvr], 'g', 'filled', 'o');
        scatter([Lhl], [Lvl], 'g', 'filled', 'v');
        scatter([Thr], [Tvr], 'b', 'filled', 'o');
        scatter([Thl], [Tvl], 'b', 'filled', 'v');
        scatter([Bhr], [Bvr], 'c', 'filled', 'o');
        scatter([Bhl], [Bvl], 'c', 'filled', 'v');
        scatter([Zhr], [Zvr], 'k', 'filled', 'o');
        scatter([Zhl], [Zvl], 'k', 'filled', 'v');
        refline(0);
        reflinexy(0,1);
        
        
        disp(num2str([Rhr / Rhl, ...
                      Lhr / Lhl, ...
                      Tvr / Tvl, ...
                      Bvr / Bvl]));  
        Trk(iN, :) = [Zhr, Zvr, Zhl, Zvl, ...
                        Rhr, Rvr, Rhl, Rvl, ...
                        Lhr, Lvr, Lhl, Lvl, ...
                        Thr, Tvr, Thl, Tvl, ...
                        Bhr, Bvr, Bhl, Bvl, ...
                        Rhr / Rhl, ...
                        Lhr / Lhl, ...
                        Tvr / Tvl, ...
                        Bvr / Bvl];
                  
        EyeTrack{iN} = {NeuronNumber, Trk(iN)};
    end
end