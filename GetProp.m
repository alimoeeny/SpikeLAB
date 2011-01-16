function [p] = GetProp(NN, FileType, StimulusType, prop)

    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NN); 
    disp(strcat('Neuron: ', num2str(NeuronNumber, '%-04.3d')));
    DataPath = GetDataPath();
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');

    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    if ( isfield(Expt.Header, prop) && ~isempty(Expt.Header.(prop)))
        p = Expt.Header.(prop);
    else if (isfield(Expt.Stimvals, prop) && ~isempty(Expt.Stimvals.(prop)))
        p = Expt.Stimvals.(prop);
        else
            p = median([Expt.Trials.(prop)]);
        end
    end
    if (length(p)>1 || isempty(p))
        disp('more than one value here');
    end