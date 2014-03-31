function [rv] = doRLSOT(iv) 

DataPath = GetDataPath();

load ../AllrlsOT.mat
AllNeurons = AllrlsOT;
clear AllrlsOT;
StimulusType = 'rls';
FileType = 'OT';

NeuronsRange = 1:length(AllNeurons)

for iN = NeuronsRange
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp([num2str(iN, '%-3.3d'), ' - Neuron: ', num2str(NeuronNumber, '%-04.3d')]);

    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    fileNames{iN} = filename;
    
    stS = expMining(Expt, [], 'st');
    ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 'bo');
    
end