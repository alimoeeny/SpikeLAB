 load('../AllTWONeurons.mat');
AllNeurons = AllTWONeurons;
clear AllTWONeurons;
FileType = 'TWO';
StimulusType = 'cylinder';



%par
for iN= [1 :length(AllNeurons)], 
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
    
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    fileNames{iN} = filename; 
    
    pD      = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
