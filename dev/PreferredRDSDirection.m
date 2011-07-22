function [pd] = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName)

    DataPath = GetDataPath();
    StimulusType = 'rds';
    FileType = 'OT';
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)) == 2)
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;
        values = unique([Expt.Trials.(Expt.Stimvals.et)]);
        [StartTime, FinishTime] = GetStartFinishTimes(FileType);
        SpikeCounts = zeros(length([Expt.Trials]),1);
        for tr = 1: length([Expt.Trials]), 
            SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
        end
        eb=[];
        for i = 1:length(values)
            eb(i) = mean(SpikeCounts([Expt.Trials(:).or]==values(i)));
            cb(i) = sum([Expt.Trials(:).or]==values(i));
        end
        pd = find(eb==max(eb));
        pd = values(pd);
    else
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.DT.mat');
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;
        pd = Expt.Stimvals.or;
    end
end
