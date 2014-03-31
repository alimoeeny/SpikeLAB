clc;
clear;
DataPath = GetDataPath();


ShowIndividualPlots = 0; % 0 or 1
SaveIndividualPlots = 0; % 0 or 1
%[AllNeurons, FileType, StimulusType] = loadAllNeurons4('DID');
[AllNeurons, FileType, StimulusType] = loadAllNeurons4('TWO');


refrenceDxs = [-0.0335   -0.0067   -0.0033         0    0.0033    0.0067    0.0335];

for iN= 1:length(AllNeurons), %iN= [20: 25] %length(AllNeurons)] %[1:20] % 1:length(AllNeurons), 
    NeuronName = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronName); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));

    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);

    TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, []);
    if(abs(TI(iN))>0.1) 
        Neuron = load(filepath);
        Expt = Neuron.Expt;
    
        values = unique([Expt.Trials(:).dx]);
        valuesLengths(iN) = length(values);
        if (length(values) == 7)
            vivals(iN, :) = values;
        else
            vivals(iN, 1) = values(1);           % min
            vivals(iN, 4) = values(ceil(end/2)); % 0
            vivals(iN, 7) = values(end);         % max
            for dxi = 2:length(values)-1
                if(values(dxi)~=0)
                    dxii = find(abs(refrenceDxs-values(dxi))==min(abs(refrenceDxs-values(dxi))));
                    disp(dxii);
                    vivals(iN, dxii) = values(dxi);
                end
            end
%         if (length(values) == 5)
%             vivals(iN, 1) = values(1); % min
%             vivals(iN, 4) = values(3); % 0
%             vivals(iN, 7) = values(5); % max
%             vivals(iN, 3) = values(2); % 0
%             vivals(iN, 5) = values(4); % 0
%         else 
%             vivals(iN, 1) = values(1);           % min
%             vivals(iN, 4) = values(ceil(end/2)); % 0
%             vivals(iN, 7) = values(end);         % max
%             vivals(iN, 3) = values(2);           % 0
%             vivals(iN, 5) = values(end - 1);     % 0
%         end
        end
    end
end