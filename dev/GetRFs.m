clear;
clc;

ShowIndividualRFs = 0;

% % % TWO
load('../AllTWONeurons.mat');
AllNeurons = AllTWONeurons;
clear AllTWONeurons;
FileType = 'TWO';
StimulusType = 'cylinder';


filenamesforbruce = {};
TI=[];


%par
for iN= 1:length(AllNeurons) %[length(AllNeurons):-1:1], 1:length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));
    switch FileType 
        case 'TWO'
            rf = RFFit(MonkeyName, NeuronNumber, ClusterName, ShowIndividualRFs);
            %TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    end
    
    %pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    RFs{iN} = {rf};
end