function [pd] = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName, DataPaths)
  
    StimulusType = 'rds';
    FileType = 'OT';
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPaths);
    if (~strcmpi(ClusterName, '.c1.'))
        if (exist(filepath, 'file')~=2)
            ClusterName = '.c1.';
            filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'OT', DataPaths);
        end
    end
   
    if(exist(filepath) == 2)
        Neuron = load(filepath);
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
        filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'DT');
        if (~strcmpi(ClusterName, '.c1.'))
            if (exist(filepath, 'file')~=2)
                ClusterName = '.c1.';
                filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'DT');
            end
        end
        if(exist(filepath) == 2)
            Neuron = load(filepath);
            Expt = Neuron.Expt;
            pd = Expt.Stimvals.or;
        else
            pd = -0009;
        end
    end
end
