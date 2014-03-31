function [SCSm, SCSse, SCSProj, conditions, varargout]= PlotSaccades(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType)
% Plots the Saccades Triggered PSTH for different conditions in the Experiment

varargout{1} = [];
varargout{2} = [];

if(nargin==1)
    Expt = MonkeyName;
    FileType = 'NNN';
    BinSize = 50;
    PlotIt = 1;
    pD = PreferredCylinderRotationDirection(Expt);
    NeuronNumber = 0;
else
    DataPath = GetDataPath();
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat');

    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;


    if strfind(FileType, 'RID')
         pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    %     if pD == -1 
    %         pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, 'DT', 0);
    %     end
    else
        pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    end
end

try
    rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName);
catch
    disp('rdsPrefDir FAILED HERE ! ! ');
    rdsPrefDir = Expt.Stimvals.or;
end

conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);

StartOffset = 500;
TrialLength = 500;


%par
for tr = 1: length([Expt.Trials]),
    avgvecX = [0];
    avgvecY = [0];
    for sc = 1: length([Expt.Trials(tr).Saccades])
        if (Expt.Trials(tr).Saccades(sc).start > StartOffset * 10) 
            %SaccadeAvgVec(tr, sc) = Expt.Stimvals.or - median(atan2(Expt.Trials(tr).Saccades(sc).pos(2), Expt.Trials(tr).Saccades(sc).pos(1)) .* (180 / (2 * pi)));
            avgvecX(end+1) = Expt.Trials(tr).Saccades(sc).pos(1);
            avgvecY(end+1) = Expt.Trials(tr).Saccades(sc).pos(2);
        end
    end
    SaccadeAvgVec(tr,1) = mean(avgvecX) ./ length(avgvecX);
    SaccadeAvgVec(tr,2) = mean(avgvecY) ./ length(avgvecY);
end
if size(SaccadeAvgVec,1)< length([Expt.Trials])
    SaccadeAvgVec(length([Expt.Trials]),1) = 0;
end


for c = 1:size(conditions,1)
    SCSm(c,:) = mean(SaccadeAvgVec(conditions(c,:),:));
    SCSse(c,:) = std(SaccadeAvgVec(conditions(c,:),:)) ./ sqrt(sum(conditions(c,:)));
    SCSProj(c,:) = sqrt(SCSm(c,1) ^ 2 + SCSm(c,2) ^ 2) .* cosd(atan2(SCSm(c,2), SCSm(c,1)) .* (90 / pi) - Expt.Stimvals.or);
    if(isnan(sum(SCSm)))
        SCSm(c,:) = [-100, -100];
        SCSse(c,:) = [-100, -100];
    end
    
end


