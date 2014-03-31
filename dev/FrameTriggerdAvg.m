
%close all, 
clear
clc

% % % % % FLID
% load('/Users/ali/DropBox/Projects/BCode/AllrdsFLIDNeurons.mat');
% AllNeurons = AllrdsFLIDNeurons;
% clear AllrdsFLIDNeurons;
% FileType = 'FLID';
% StimulusType = 'cylinder';

% % % % DRID
% load('/Users/ali/DropBox/Projects/BCode/AllrdsDRIDNeurons.mat');
% AllNeurons = AllrdsDRIDNeurons;
% clear AllrdsDRIDNeurons;
% FileType = 'dRID';
% StimulusType = 'rds';

% % % % rdsDT
load('../AllrdsDTNeurons.mat');
AllNeurons = AllrdsDTNeurons;
clear AllrdsDTNeurons;
FileType = 'DT';
StimulusType = 'rds';

% % % cylinderDT
% load('../AllcylinderDTNeurons.mat');
% AllNeurons = AllcylinderDTNeurons;
% clear AllcylinderDTNeurons;
% FileType = 'DT';
% StimulusType = 'cylinder';

% % DID
% load('../AllDIDNeurons.mat');
% AllNeurons = AllDIDNeurons;
% clear AllDIDNeurons;
% FileType = 'DID';
% StimulusType = 'cylinder';
 

% % % % rdsDT IN DID CELLS
% load('/Users/ali/DropBox/Projects/BCode/AllDIDrdsDTNeurons.mat');
% % %%load('/Users/ali/DropBox/Projects/BCode/AllrdsDTNeurons.mat');
% AllNeurons = AllDIDrdsDTNeurons;
% clear AllDIDrdsDTNeurons;
% FileType = 'DT';
% StimulusType = 'rds';


% % % % cylinderDT IN DID CELLS
% load('/Users/ali/DropBox/Projects/BCode/AllDIDrdsDTNeurons.mat');
% % % %%load('/Users/ali/DropBox/Projects/BCode/AllrdsDTNeurons.mat');
% AllNeurons = AllDIDrdsDTNeurons;
% clear AllDIDrdsDTNeurons;
% FileType = 'DT';
% StimulusType = 'rds'; %'cylinder';
% 

kernel = gausswin(9);
frameTrig = [];
for iN = 1:length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[c, eb] = PlotISISdf(MonkeyName, NeuronNumber, ClusterName, FileType, 'cylinder', 1);
    %[c, eb] = PlotISISdf(MonkeyName, NeuronNumber, ClusterName, FileType, 'rds', 1);
    %continue 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DataPath = GetDataPath();
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)))
    %if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName,      'rds.' , FileType,'.mat'))) && ... 
    %    exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'cylinder.' , FileType,'.mat'))))
        
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;
        FrameDur = 0;
        if (isfield(Expt.Header, 'frameperiod')) 
            FrameDur = Expt.Header.frameperiod/10;
        end
        if (FrameDur < 5)
            disp(['FrameDur problem! ', num2str(FrameDur)]);
            FrameDur = 16.70100;
        end
        if(isfield(Expt.Trials,'dur'))
            triallength = round(median([Expt.Trials.dur]) / 10);
        else
            triallength = round((Expt.Trials(1).End - Expt.Trials(1).Start)/10);
        end
        frames = triallength/FrameDur;
        trialCount = length(Expt.Trials);
        cntr = 0; acr = [];
        Spks = zeros(round(frames*trialCount), triallength);
        Stimuli = zeros(round(frames*trialCount),1);
        for trialNumber = 1:trialCount
            spikeTimes = round(0.1*[Expt.Trials(trialNumber).Spikes]);
            spikeTimes = spikeTimes(spikeTimes>0 & spikeTimes < triallength);
            for fn = 1:frames
                %spks = spikeTimes(spikeTimes>= (fn-1)*FrameDur & spikeTimes< fn*FrameDur);
                spks = round(spikeTimes - ((fn-1)*FrameDur));
                spks = spks(spks > 0);
                if(sum(spks)>0)
                    cntr = cntr +1;
                    Spks(cntr, spks) = 1;
                    Stimuli(cntr) = Expt.Trials(trialNumber).(Expt.Stimvals.et);
                    if ((Stimuli(cntr) == 0) && strcmpi(FileType, 'DID'))
                        Stimuli(cntr) = 10 * Expt.Trials(trialNumber).Id;
                    end
                end
            end
            %disp(cntr);
        end
        frameTrig = [];
        spikes = convn(Spks, kernel');
        values = unique(Stimuli);
        for s = 1:length(values)
            frameTrig(s,:) = mean(spikes(Stimuli==values(s),:));
        end
        pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
        TI = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
        SF{iN} = {frameTrig, values, pD, TI};
    end
end

%%
fTPop = [];
for iN = 1: length(SF)
    frameTrig = SF{iN}{1};
    values = SF{iN}{2};
    pD = SF{iN}{3};
    switch FileType
        case {'DID'}
            start=1; finish = 1900;
        case {'DT'}
            start=1; finish = 350;
    end
    if (pD == 2)
        fTPop(iN, :, :) = [frameTrig(1,start:finish); frameTrig(2,start:finish); frameTrig(3,start:finish); frameTrig(end-2,start:finish); frameTrig(end-1,start:finish); frameTrig(end,start:finish)];
    else
        fTPop(iN, :, :) = [frameTrig(end,start:finish); frameTrig(end-1,start:finish); frameTrig(end-2,start:finish); frameTrig(3,start:finish); frameTrig(2,start:finish); frameTrig(1,start:finish)];
    end
end

%% Graphics

figure(1222), clf, hold on
h = plot(squeeze(mean(fTPop,1))');
xlim([0 250]);
set(h, 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xtl = [0, 50, 100, 150, 200, 250];
set(gca, 'XTick', xtl);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);

%%
for iN = 1: length(SF)
    figure(1118+iN)
    frameTrig = SF{iN}{1};
    values = SF{iN}{2};
    h = plot(frameTrig');
    xlim([0 250]);
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
    xtl = [0, 50, 100, 150, 200, 250];
    set(gca, 'XTick', xtl);
    set(gca, 'XTickLabel', {num2str(xtl')});
    legend(h, num2str(values));
    title(AllNeurons(iN));
    pause;
end

figure(1118), clf, hold on,
plot(frameTrig');

