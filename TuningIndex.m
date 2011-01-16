function [TI] = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType, reqparam)
% Gets the Tuning Index of the cell  based on the average firing rate  

%Prep

DataPath = GetDataPath();

FileType = ExperimentType;
switch FileType
    case {'SRID', 'DRID'}
        %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.','orXId.', FileType,'.mat');
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.','DT.mat');
        if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~= 2)
            filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'XIdD.DT.mat');
        end
    case 'BDID'
%         filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.DT.mat');
%         if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~= 2)
%             filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.DIDB.mat');
%             if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~=2)
%                 filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.DID.mat');
%             end
%         end
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');

    otherwise
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    
end

Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
Expt = Neuron.Expt;

if(nargin<6)
    if strcmp(FileType,'BDID')
         reqparam = 'Id'; %reqparam = 'dx';
    else
        reqparam = Expt.Stimvals.et;
    end
end
    
switch FileType
    case 'TWO'
        StartTime = 1000;
        values = unique([Expt.Trials.('dx')]);
    case 'DID'
        StartTime = 500;
        values = unique([Expt.Trials.(reqparam)]);
    case 'BDID'
        StartTime = 500;
        values = unique([Expt.Trials.(reqparam)]);
    otherwise
        StartTime = 500;
        values = unique([Expt.Trials.(reqparam)]);
end

if(isfield(Expt.Trials, 'dur'))
    FinishTime = round(median([Expt.Trials(:).dur])) + 500;
else
    FinishTime = round(median([Expt.Trials(:).End] - [Expt.Trials(:).Start])) + 500;
end
if FinishTime<6000
    StartTime = 500;
end
if (strcmp(FileType, 'BDID')|strcmp(FileType, 'DID'))
    StartTime = 500;
    FinishTime = 5500;
end
SpikeCounts = zeros(length([Expt.Trials]),1);
for tr = 1: length([Expt.Trials]), 
    SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
end

eb=[];
for i = 1:length(values)
    if strcmp(FileType,'TWO')
        eb(i) = mean(SpikeCounts([Expt.Trials.('dx')]==values(i)));
        cb(i) = sum([Expt.Trials.('dx')]==values(i));
    else
        eb(i) = mean(SpikeCounts([Expt.Trials.(reqparam)]==values(i)));
        cb(i) = sum([Expt.Trials.(reqparam)]==values(i));
    end
end

TI = (sum(eb(values>0) .* cb(values>0)) ./ mean(cb(values>0)) - sum(eb(values<0) .* cb(values<0)) ./ mean(cb(values<0)) ) / ...
     (sum(eb(values>0) .* cb(values>0)) ./ mean(cb(values>0)) + sum(eb(values<0) .* cb(values<0)) ./ mean(cb(values<0)) );










