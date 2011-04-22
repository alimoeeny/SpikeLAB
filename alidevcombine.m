function r = alidevcombine(NeuronName)
    DataPath = GetDataPath();
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus({NeuronName});
    filename = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), '.mat');
    cd /bgc/bgc/matlab/dev/
    combine(filename);
    r = 0;