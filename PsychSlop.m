function [slp, hf] = PsychSlop(NeuronName, StimulusType, ExperimentType, reqparam, showthefit)
  
if (nargin>3)
    %Prep
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronName);
    DataPath = GetDataPath();

    FileType = ExperimentType;

%    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);

    Neuron = load(filepath);
    Expt = Neuron.Expt;

% if (isempty(reqparam))
%     reqparam = Expt.Stimvals.et;
    if strcmpi(reqparam, 'TwoCylDisp')
        reqparam = 'dx';
    end
% end

else
    Expt = NeuronName;
    if ~(strcmp(StimulusType,''))
        reqparam = StimulusType;
        StimulusType = '';
    else
        reqparam = 'dx';
    end
    if (strcmpi(ExperimentType,''))
        ExperimentType = 'DID';
    end
end


values = unique([Expt.Trials.(reqparam)]);
if strcmp(reqparam, 'dx')
    values = values(values<0.5 & values>-0.5);
    tmpdxvals = [];
    for dxx = 1: length(values)
       if (sum([Expt.Trials(:).dx]==values(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4 && sum([Expt.Trials(:).dx]== -values(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4)
           tmpdxvals(end+1) = values(dxx);
       end
    end
    values = tmpdxvals;     
end

switch ExperimentType
    case {'DID', 'ABD', 'TWO', 'DIDB'}
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
    case 'BDID'
        if(mean([Expt.Trials([Expt.Trials(:).bd]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
end

Responses = []; %zeros(length(values),1);
for tr = 1: length(values), 
    Responses(tr).x = values(tr);
%    Responses(tr).r    = sum([Expt.Trials([Expt.Trials(:).dx]==values(tr)).RespDir]==ResponseToPositive) / sum([Expt.Trials([Expt.Trials(:).dx]==values(tr)).RespDir]~=0);
    Responses(tr).resp = sum([Expt.Trials([Expt.Trials(:).(reqparam)]==values(tr)).RespDir]==ResponseToPositive);
    Responses(tr).n =  sum([Expt.Trials([Expt.Trials(:).(reqparam)]==values(tr)).RespDir]~=0);
end
         
hf = 0;
if showthefit
    hf = figure(18186); clf, hold on,
    psf = fitpsf(Responses, 'showfit');
else
    psf = fitpsf(Responses);
end
tempResponses = [];
for i = 1:length(psf.data)
    if (psf.data(i).p<0.9 & psf.data(i).p>0.1) 
        sk = length(tempResponses)+1;
        tempResponses(sk).x = Responses(i).x;
        tempResponses(sk).resp = Responses(i).resp;
        tempResponses(sk).n = Responses(i).n;
    end
end
if (length(tempResponses)>2)
   if showthefit
       psf = fitpsf(tempResponses, 'showfit');
   else
    psf = fitpsf(tempResponses);
   end
end
slp = psf;

