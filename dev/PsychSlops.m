function [slp] = PsychSlops(NeuronNumber, ClusterName, StimulusType, ExperimentType, reqparam)
  
if(nargin==1)
    Expt = NeuronNumber;
else if nargin>2
        %Prep
        MonkeyName = 'icarus';
        DataPath = '/bgc/data/';

        FileType = ExperimentType;

        filename = strcat('ic', num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');

        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;
    end
end

if(nargin<5)
    reqparam = Expt.Stimvals.et;
end

for blks = 1: length(Expt.Header.BlockStart)
    if length(Expt.Header.BlockStart) == 1
        trls = 1:length(Expt.Trials);
    else
        if blks < length(Expt.Header.BlockStart)
            trls = (([Expt.Trials(:).Trial] >= Expt.Header.BlockStart(blks)) & [Expt.Trials(:).Trial] < Expt.Header.BlockStart(blks+1));
        else
            trls = (([Expt.Trials(:).Trial] >= Expt.Header.BlockStart(blks)) & [Expt.Trials(:).Trial] < Expt.Trials(end).Trial );
        end
    end
    values = unique([Expt.Trials(trls).(reqparam)]);

    if (sum(trls) > length(values) * 4)
        if strcmp(reqparam, 'bd')
            if(mean([Expt.Trials([Expt.Trials(:).bd]>0).RespDir])>0) % no need for trls
                    ResponseToPositive = 1;
                    ResponseToNegative = -1;
            else
                    ResponseToPositive = -1;
                    ResponseToNegative = 1;
            end
        
        else
            if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0) % no need for trls
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
        %    Responses(tr).r    = sum([Expt.Trials([Expt.Trials(trls).dx]==values(tr)).RespDir]==ResponseToPositive) / sum([Expt.Trials([Expt.Trials(trls).dx]==values(tr)).RespDir]~=0);
            Responses(tr).resp = sum([Expt.Trials([Expt.Trials(trls).(reqparam)]==values(tr)).RespDir]==ResponseToPositive);
            Responses(tr).n =  sum([Expt.Trials([Expt.Trials(trls).(reqparam)]==values(tr)).RespDir]~=0);
        end

        tempResponses = [];
        if (min([Responses.n])>0)
            %psf = fitpsf(Responses, 'showfit');
            psf = fitpsf(Responses);
            if(psf.fit(2)<100)
                for i = 1:length(Responses)
                    if (psf.data(i).p<0.9 & psf.data(i).p>0.1) 
                        sk = length(tempResponses)+1;
                        tempResponses(sk).x = Responses(i).x;
                        tempResponses(sk).resp = Responses(i).resp;
                        tempResponses(sk).n = Responses(i).n;
                    end
                end
            else
                pst = [];
            end
            if (length(tempResponses)>2)
                psf = fitpsf(tempResponses);
            end
        
        slp(blks) = psf;
        end
    end
end
