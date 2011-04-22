function [p] = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, varargin)

    DataPath = GetDataPath();

    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    
    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    
    or = Expt.Stimvals.or;
    if isempty(or)
        or = median([Expt.Trials.or]);
    end
    p = [or];