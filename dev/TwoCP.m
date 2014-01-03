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
    
    
    PrefCyldx = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);

    if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
        ResponseToPositive = 1;
        ResponseToNegative = -1;
    else
        ResponseToPositive = -1;
        ResponseToNegative = 1;
    end
    
    COND = logical([]);
    if(PrefCyldx == 2)
        COND(1,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
        COND(2,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
        
        COND(3,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
        COND(4,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
    else
        COND(1,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
        COND(2,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;

        COND(3,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
        COND(4,:) = [Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
    end
    
    
    SpikeCounts = zeros(length([Expt.Trials]),1);
    for tr = 1: length([Expt.Trials]),
        SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
    end

    
    CP(iN) = ROCAUC(SpikeCounts(COND(1,:)), SpikeCounts(COND(2,:)));
    CPp(iN)= ROCAUC(SpikeCounts(COND(3,:)), SpikeCounts(COND(4,:)));
        
end

figure(1134), clf, hold on
plot(CP, 'b');
plot(CPp, 'r');
refline(0,0.5);

figure(2455), clf, hold on
scatter(CP, CPp, 'filled');
refline(1);

