%

clear;

[AllNeurons, FileType, StimulusType] = loadAllNeurons4('TWO');

for iN= 1 : length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);

    disp(filename);
    
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    
    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName);
    conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);

    StartTime =  500;
    FinishTime = 20000;

    SpikeRates = zeros(length([Expt.Trials]),1);
    for tr = 1: length([Expt.Trials]),
        SpikeRates(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
        SpikeRates(tr) = SpikeRates(tr) ./ ((FinishTime-StartTime)/10000);
    end
    
    SpikeRates = zscore(SpikeRates);

    values = unique([Expt.Trials.('dx')]);
     
    eb = []; cb = [];
    for cnd = 1: size(conditions, 1)
        eb(cnd) = mean(SpikeRates(conditions(cnd,:)));
        cb(cnd) = std(SpikeRates(conditions(cnd,:)))./sqrt(sum(conditions(cnd,:)));
    end
    
    for i = 1:length(values)
        eb(end+1) = mean(SpikeRates([Expt.Trials(:).dx]==values(i)));
        cb(end+1) = std(SpikeRates([Expt.Trials(:).dx]==values(i)))./sqrt(sum([Expt.Trials(:).dx]==values(i)));
    end

    
    %% Load the cylinder.DT
   try 
    filepathDT = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'Cylinder', 'DT');
    
    NeuronDT = load(filepathDT);
    ExptDT = NeuronDT.Expt;
    
    StartTimeDT =  500;
    FinishTimeDT = 5000;

    SpikeRatesDT = zeros(length([ExptDT.Trials]),1);
    for tr = 1: length([ExptDT.Trials]),
        SpikeRatesDT(tr) = sum([ExptDT.Trials(tr).Spikes]>=StartTimeDT & [ExptDT.Trials(tr).Spikes]<=FinishTimeDT);
        SpikeRatesDT(tr) = SpikeRatesDT(tr) ./ ((FinishTimeDT-StartTimeDT)/10000);
    end
    
    SpikeRatesDT = zscore(SpikeRatesDT);
    
    valuesDT = unique([ExptDT.Trials.('dx')]);
    for i = 1:length(valuesDT)
        eb(end+1) = mean(SpikeRatesDT([ExptDT.Trials(:).dx]==valuesDT(i)));
        cb(end+1) = std(SpikeRatesDT([ExptDT.Trials(:).dx]==valuesDT(i)))./sqrt(sum([ExptDT.Trials(:).dx]==valuesDT(i)));
    end
    

    
    
    
    %%
    h = figure(7788); clf; hold on;
        
%   errorbar(valuesDT, eb(end-length(valuesDT)+1:end), cb(end-length(valuesDT)+1:end), 'k');

%    errorbar(values, eb([size(conditions,1)+1:size(conditions,1) + length(values)]), cb([size(conditions,1)+1:size(conditions,1)+length(values)]))
    
    mmm = mean(eb(size(conditions,1) + length(values) + find(valuesDT<0)));
    mmmc = std(eb(size(conditions,1) + length(values) + find(valuesDT<0))) ./ sqrt(sum(size(conditions,1) + length(values) + find(valuesDT<0)));
    mmp = mean(eb(size(conditions,1) + length(values) + find(valuesDT>0)));
    mmpc = std(eb(size(conditions,1) + length(values) + find(valuesDT>0))) ./ sqrt(sum(eb(size(conditions,1) + length(values) + find(valuesDT>0))));
    
    x = [min(values) -0.0 0.0 max(values)];
    y = [mmm eb(3), eb(6) mmp];
    e = [mmmc cb(3), cb(6) mmpc];

    errorbar(x,y,e, 'r', 'LineWidth', 4, 'LineStyle', 'none')
    
    
    
    
    %%
   
    R(iN,:) = [y, e];
    
   catch e
       disp(e.message)
   end
    
end
