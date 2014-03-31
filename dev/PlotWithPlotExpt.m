clear;
clc
cd /Users/ali/Dropbox/Projects/SpikeLAB/dev
DataPath = GetDataPath();

ShowSingleCellSDFs = 0; % 0 or 1
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('TWO');

for iN= [1:length(AllNeurons)] %[1:33 40:length(AllNeurons)],
    
    IdColor{iN} = [1 0.5 0.1];
    DotSizes(iN) = 100;
    toolowFR = 0;
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber);
    TI(iN)      = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, []);
    TIcyldx(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, 'DT',     []);
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
    
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    %    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    fileNames{iN} = filename;
    
    figure(1112), clf, 
    hp = PlotExpt(Expt);
    
    print(hp.fig, '-dpsc', '-r300', '-zbuffer', ['../figs/PSTH', '-', date, '-', filename, '.eps']);
    
end