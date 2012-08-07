function [p] = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, varargin)

    DataPath = GetDataPath();

%    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    
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
        end
    end