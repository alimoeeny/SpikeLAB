function fpath = MakeFilePath (MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType)
    DataPaths = GetDataPath(); 
    if ClusterName(2) == 'e'
         fname = [MonkeyAb(MonkeyName), 'M', num2str(NeuronNumber, '%-04.3d'), '.', StimulusType,'.', ExperimentType, '.', ClusterName, '.mat'];
         fpath = [DataPaths{2}, MonkeyName, '/M', num2str(NeuronNumber, '%-04.3d'), '/' ,fname];
     else
         fname = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', ExperimentType,'.mat');
         fpath = [DataPaths{1}, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , fname];
     end