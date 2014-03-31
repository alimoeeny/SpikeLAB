function r = alicombine(NeuronName)
    DataPath = '/sc/bgc/bgc/data/';
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus({NeuronName});
    filename = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), '.mat');
    combine(filename);
    r = 0;