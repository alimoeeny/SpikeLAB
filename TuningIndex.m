function [TI] = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType, reqparam, doitsquare)
% Gets the Tuning Index of the cell  based on the average firing rate  

if doitsquare
    disp('going Square');
end

DataPaths = GetDataPath();

FileType = ExperimentType;
switch FileType
    case {'SRID', 'DRID'}
        %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.','orXId.', FileType,'.mat');
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.','DT.mat');
        if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~= 2)
            filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'XIdD.DT.mat');
        end

    otherwise
        filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType);
        filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType);
end

if ((exist(strcat(DataPaths{1}, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~= 2) ...
    && (exist(strcat(DataPaths{2}, MonkeyName, '/M', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) ~= 2))
    TI = -999;
    return
end

Neuron = load(filepath);
Expt = Neuron.Expt;

if( (nargin<6) || isempty(reqparam))
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
%        values = values(abs(values) > 0.01);
%        values = values(abs(values)< 0.01);
%         a = floor(length(values)/2);
%         if(a < length(values) / 2)
%             b = a + 2;
%         else
%             b = a + 1;
%         end
%         values = values([a b]);
%        disp('TI BASED ONLY ON VALUES BETWEEN -0.01 and 0.01 - - - - -- - - - - - - ');
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
if (strcmp(FileType, 'BDID')||strcmp(FileType, 'DID'))
    StartTime = 500;
    FinishTime = 5500;
end
SpikeCounts = zeros(length([Expt.Trials]),1);
for tr = 1: length([Expt.Trials]), 
    if doitsquare
        SpikeCounts(tr) = sqrt(sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime));
    else
        SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
    end
end

eb=[]; cb = [];
for i = 1:length(values)
    if strcmp(FileType,'TWO')
        eb(i) = mean(SpikeCounts([Expt.Trials.('dx')]==values(i)));
        cb(i) = sum([Expt.Trials.('dx')]==values(i));
    else
        eb(i) = mean(SpikeCounts([Expt.Trials.(reqparam)]==values(i)));
        cb(i) = sum([Expt.Trials.(reqparam)]==values(i));
    end
end

tis = [];
for i = 1:floor(length(eb)/2)
    a = i;
    b = 1 + length(eb) - i;
    tis(i) = (sum(eb(a) .* cb(a)) ./ mean(cb(a)) - sum(eb(b) .* cb(b)) ./ mean(cb(b)) ) / ...
     (sum(eb(a) .* cb(a)) ./ mean(cb(a)) + sum(eb(b) .* cb(b)) ./ mean(cb(b)) );
end

TI = max(tis);
if (TI < abs(min(tis)))
    TI = min(tis);
end

if(isempty(TI))
    TI = -998;
end

% TI = (sum(eb(values>0) .* cb(values>0)) ./ mean(cb(values>0)) - sum(eb(values<0) .* cb(values<0)) ./ mean(cb(values<0)) ) / ...
%      (sum(eb(values>0) .* cb(values>0)) ./ mean(cb(values>0)) + sum(eb(values<0) .* cb(values<0)) ./ mean(cb(values<0)) );









