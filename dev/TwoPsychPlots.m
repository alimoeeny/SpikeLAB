clear;
clc;

[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('TWO');
%Prep
DataPaths = GetDataPath();


for iN= 1 : length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN));
    if(isempty(ClusterName)), ClusterName = ['cell' num2str(AllNeurons{iN,2})]; end
    disp(['iN: ' ,num2str(iN) , ' , Neuron: ', MonkeyName, ' - ',  num2str(NeuronNumber, '%-04.3d'), ' - ', ClusterName]);
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    clear Neuron
    
    clf, hold on,
    ExptPsych(Expt);
    hp = gcf();
    print(hp, '-dpsc', '-r300', '-zbuffer', ['../figs/TWOPsych', '-', date, '-', filename, '.eps']);
        
end


