function [p] = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, varargin)

if(nargin>2)    
    DataPath = GetDataPath('server');
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);
    
    Neuron = load(filepath);
    Expt = Neuron.Expt;
else
    Expt = MonkeyName;
end
    try
        disp(unique([Expt.Trials.or]));
        disp('--------------------------------------------------------');
    catch e
        
    end
    or = Expt.Stimvals.or;
    if isempty(or)
        or = median([Expt.Trials.or]);
    end
    p = [or];
    
    if(~isempty(varargin))
        if(strcmp(varargin{1},'ve'))
            p = Expt.Stimvals.ve;
        else 
            p = Expt.Stimvals.(varargin{1});
            if isempty(p)
                p = median([Expt.Trials.(varargin{1})]);
            end
        end
    end