clear, clc

DataPath = GetDataPath('server');

% % % % DID
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('DID');

for iN = 1:length(AllNeurons)
   [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN));
   disp([num2str(iN), ' - ', num2str(NeuronNumber)]);
   
   filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);
   filepath = [filepath(1:end-3), 'EYE.mat'];
   
   if(exist(filepath,'file'))
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    
    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, DataPath, FileType, 0);
    try
        rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName, DataPath);
    catch
        disp('rdsPrefDir FAILED HERE ! ! ');
        rdsPrefDir = Expt.Stimvals.or;
    end
    conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);
    
    StartTime = 1;
    FinishTime = 1500;
    
    
    EyeTraces = [];
    for i = 1: length(Expt.Trials)
        EyeTraces(i, :, :) = Expt.Trials(i).EyeData;
    end
    
    %AvgTracks(iN, :, :) = squeeze(mean(EyeTraces));
    
    a = squeeze(mean(EyeTraces(:,:,:)));
    HvergenceBaseline = a(:,1) - a(:,3);
    VvergenceBaseline = a(:,1) - a(:,3);
    
    a1 = squeeze(mean(EyeTraces(conditions(1,:), :, :)));
    Hvergence1 = a1(:,1) - a1(:,3);
    Vvergence1 = a1(:,2) - a1(:,4);
    
    a2 = squeeze(mean(EyeTraces(conditions(2,:), :, :)));
    Hvergence2 = a2(:,1) - a2(:,3);
    Vvergence2 = a2(:,2) - a2(:,4);
    
    if(size(Hvergence1,1)>=FinishTime)
        ChanCorrs(iN) = corr(a(1:1200,1), a(1:1200,2));
        ChanVars(iN, :) = [var(a(1:1200,1)), var(a(1:1200,2))];
        ViVigance(iN, :,:) =  [Hvergence1(StartTime:FinishTime), Hvergence2(StartTime:FinishTime), HvergenceBaseline(StartTime:FinishTime)];
    else
        disp (size(Hvergence1));
    end
    
   else
    disp(['file not found', filepath]);   
   end
end


%%
for i = 1:size(ViVigance,2)
    hV1(i) = mean(ViVigance(:,i,1) - ViVigance(:,i,3));
    hV2(i) = mean(ViVigance(:,i,2) - ViVigance(:,i,3));
    th(i) = ttest(ViVigance(:,i,1) - ViVigance(:,i,2));
end

figure, plot(hV1)
hold on, plot(hV2)
