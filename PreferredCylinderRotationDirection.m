function [pd] = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, varargin)
% 1: positive, 2: negative , -1 when file does not exist
% arguments Monkey name,  neuron number, clustername, fiiletype, Include zero or not,
% 'forced' (forced to use the filetype given not to use ABD if available or ...)

% manual mapping for tricky neurons
% positive pref
if (sum(strcmpi([MonkeyName, num2str(NeuronNumber)],{'dae523'; ''}))>0)
    pd = 1;
    return
end
% negative pref
if (sum(strcmpi([MonkeyName, num2str(NeuronNumber)],{''; ''}))>0)
    pd = 2;
    return
end


if(nargin==1)
    Expt = MonkeyName;
    FileType = 'NNN';
    IncludeZero = 0;
else

    DataPath = GetDataPath();
    StimulusType = 'cylinder';

    if(size(varargin,2)>0)
        FileType = varargin{1};
    else
        FileType = 'ABD';
    end
    

    if (size(varargin,2) > 2)
        if(strcmpi(varargin{3}, 'forced')==0)
           if strcmpi(FileType, 'DPI')
                filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', 'ABD.mat');
                if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename),'file')==2)
                    FileType = 'ABD';
                end
           end
        end
    end
    if ~isempty(strfind(FileType, 'RID'))
        StimulusType = 'rds';
    end

    if(size(varargin,2)>1)
        IncludeZero = varargin{2};
    else
        IncludeZero = 1;
    end

    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    
%     if ~isempty(strfind(FileType, 'RID'))
%         filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', 'DRID','.mat');
%         if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename),'file')~=2)
%             filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'cylinder.ABD','.mat');
%         end
%     end

    if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename),'file')~=2)
        pd = -1;
        disp(['FILE NOT FOUND - THERE IS A PROBLEM HERE! ! !', strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)]);
        return;
    end

    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
end


if strcmp(FileType,'BDID') 
    values = unique([Expt.Trials.('Id')]);
else if strcmp(FileType,'TWO') 
    values = unique([Expt.Trials.('dx')]);
else if ~isempty(strfind(FileType, 'RID'))
        values = unique([Expt.Trials(:).Id]);
    else
        values = unique([Expt.Trials.(Expt.Stimvals.et)]);
    end
    end
end

if strcmp(FileType, 'DID') || strcmp(FileType, 'BDID')
    StartTime = 500;
    FinishTime = 5500;
else
    [StartTime, FinishTime] = GetStartFinishTimes(FileType);
end

SpikeCounts = zeros(length([Expt.Trials]),1);
for tr = 1: length([Expt.Trials]), 
    SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
end

eb=[];
for i = 1:length(values)
    if strcmpi(FileType, 'BDID')
        eb(i) = mean(SpikeCounts([Expt.Trials(:).Id]==values(i)));
        cb(i) = sum([Expt.Trials(:).Id]==values(i));
    else
        if ~isempty(strfind(FileType, 'RID'))
            eb(i) = mean(SpikeCounts([Expt.Trials(:).Id]==values(i)));
            cb(i) = sum([Expt.Trials(:).Id]==values(i));
        else
            eb(i) = mean(SpikeCounts([Expt.Trials(:).dx]==values(i)));
            cb(i) = sum([Expt.Trials(:).dx]==values(i));
        end
    end
end

pos(1) = sum(eb(values>0) .* cb(values>0)) ./ sum(cb(values>0));
pos(2) = sum(eb(values<0) .* cb(values<0)) ./ sum(cb(values<0));
if(IncludeZero)
    pos(3) = sum(eb(values==0) .* cb(values==0)) ./ sum(cb(values==0));
end

if pos(1) == pos(2)
    pd = -1;
    disp('We are in a pickle!  ********************************************************************');
else
    pd = find(pos==max(pos));
end
