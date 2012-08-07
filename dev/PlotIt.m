[AllNeurons, FileType, StimulusType] = loadAllNeurons4('TWO');

for iN= 1 : length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);

    disp(filename);
    
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    
    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName);
    conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);

    
end
