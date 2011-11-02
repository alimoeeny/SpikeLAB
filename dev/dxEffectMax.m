function [dxMax] = dxEffectMax(MonkeyName, NeuronNumber, ClusterName, StimulusType)

DataPath = GetDataPath();

filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.DT.mat');
if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', filename)) ~= 2)
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), '.c1.', StimulusType,'.DT.mat');
    if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', filename)) ~= 2)
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName,'cylinder.DT.mat');
    end
end
Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
Expt = Neuron.Expt;

StartTime = 500;
FinishTime = median([Expt.Trials(:).End] - [Expt.Trials(:).Start])+StartTime;

dxValues = unique([Expt.Trials.dx]);

SpikeCounts = zeros(length([Expt.Trials]),1);
for tr = 1: length([Expt.Trials]), 
    SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
end


eb=[]; cb = [];
for i = 1:length(dxValues)
    eb(i) = mean(SpikeCounts([Expt.Trials(:).dx]==dxValues(i)));
    cb(i) = sum([Expt.Trials(:).dx]==dxValues(i));
end

dxMax = (max(eb) - min(eb)) ./ ((FinishTime - StartTime) ./ 10000);


end
