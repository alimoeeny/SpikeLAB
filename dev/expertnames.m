for i = 1:length(AllPursuitNeurons), 
    DataPath = GetDataPath();
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus({AllPursuitNeurons{i}});
    filename = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), '.mat');
    disp(filename);
end