clear;
clc
cd /Users/ali/Dropbox/Projects/SpikeLAB/dev
DataPath = GetDataPath();

ShowSingleCellSDFs = 0; % 0 or 1
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('TWO');



for iN= [1:length(AllNeurons)]
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber);
    disp([num2str(iN), ' - ' , MonkeyName, num2str(NeuronNumber), ClusterName]);
    
    TI(iN)      = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, []);
    
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    
    Neuron = load(filepath); Expt = Neuron.Expt; clear Neuron;
    fileNames{iN} = filename;
    
    if(TI(iN)>0)
        pD = 1;
    elseif TI(iN)<0
        pD = 2;
    else
        pD = 0;
    end
    
    conditions = logical([]);
    if strcmp(FileType, 'DID')
        conditions = GetConditions(Expt, FileType, pD);
    elseif strcmp(FileType, 'TWO')
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        if(pD == 2)
            conditions(1,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(4,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(5,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(6,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(7,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(8,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare
            conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(10,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(11,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            
            % 16 - 19 are the the two dxs near zero disparity.
            conditions(16,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(17,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(18,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(19,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            
            % 24 - 27 are the the two dx zero disparitis for Z-Scored CP calculation.
            conditions(24,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(25,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(26,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(27,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            
            
        else
            conditions(1,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(2,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
            conditions(3,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(4,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]~=0;
            conditions(5,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(6,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(7,:) = [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(8,:) = [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip pair to copmare
            conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]~=0;
            conditions(10,:)= [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip
            conditions(11,:)= [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; %flip pair to compare
            
            % 16 - 19 are the the two dxs near zero disparity.
            conditions(16,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(17,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(18,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(19,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
            
            
            % 24 - 27 are the the two dx zero disparitis for Z-Scored CP calculation.
            conditions(24,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(25,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            conditions(26,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
            conditions(27,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
            
        end
        
    end
    stepSize =  500;
    winWidth = 1000;
    
    clear GCP
    parfor iWin = 0:20000 / stepSize
        
        StartTime = iWin * stepSize;
        FinishTime = StartTime + winWidth;
        
        SpikeCounts = zeros(length([Expt.Trials]),1);
        for tr = 1: length([Expt.Trials]),
            SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
        end
        
        
        % FINAL CP
        zpId = zscore([SpikeCounts(conditions(24,:)); SpikeCounts(conditions(25,:))]);
        znId = zscore([SpikeCounts(conditions(26,:)); SpikeCounts(conditions(27,:))]);
        pdx  = zscore([SpikeCounts(conditions(16,:)); SpikeCounts(conditions(17,:))]);
        ndx  = zscore([SpikeCounts(conditions(18,:)); SpikeCounts(conditions(19,:))]);
        
        mtn = 6; % minimum trials needed
        aa = []; bb = []; AA={}; BB={};
        % basic CP in zero disparity trials with null id
        if ((sum(conditions(24,:))>=mtn) && ((sum(conditions(25,:))>=mtn) )) ,
            at = zpId(1:sum(conditions(24,:)));
            bt = zpId(sum(conditions(24,:))+1:sum(conditions(24,:))+sum(conditions(25,:)));
            aa = [aa; at];
            bb = [bb; bt];
            AA{end+1} = at;
            BB{end+1} = bt;
        end
        % basic CP in zero disparity trials with pref id
        if ((sum(conditions(26,:))>=mtn) && ((sum(conditions(27,:))>=mtn) )) ,
            at = znId(1:sum(conditions(26,:)));
            bt = znId(sum(conditions(26,:))+1:sum(conditions(26,:))+sum(conditions(27,:)));
            aa = [aa; at];
            bb = [bb; bt];
            AA{end+1} = at;
            BB{end+1} = bt;
        end
        %CP in next to zero dx trials with null DX (no id here)
        if ((sum(conditions(16,:))>=mtn) && ((sum(conditions(17,:))>=mtn) )) ,
            at = pdx(1:sum(conditions(16,:)));
            bt = pdx(sum(conditions(16,:))+1:sum(conditions(16,:))+sum(conditions(17,:)));
            aa = [aa; at];
            bb = [bb; bt];
            AA{end+1} = at;
            BB{end+1} = bt;
        end
        %CP in next to zero dx trials with pref DX (no id here)
        if ((sum(conditions(18,:))>=mtn) && ((sum(conditions(19,:))>=mtn) )) ,
            at = ndx(1:sum(conditions(18,:)));
            bt = ndx(sum(conditions(18,:))+1:sum(conditions(18,:))+sum(conditions(19,:)));
            aa = [aa; at];
            bb = [bb; bt];
            AA{end+1} = at;
            BB{end+1} = bt;
        end
        
        GCP(iWin+1) = ROCAUC(aa, bb);
        [roc2compare, p, pp] = ROCAUCSignificanceGrandCP(AA, BB);
        %GCPSig(iWin+1) = p;
        
    end
    GrandCP(iN, :) = GCP;
    %GrandCPSig(iN,:) = GCPSig;
end

%% save it
save(['~/Desktop/matlab.SlidingGrandCP.', FileType, num2str(now) '.mat'], '-v7.3')


%% average it all
for j = 1:size(GrandCP, 2)
    ma = [];
    for i = 1:size(GrandCP, 1)
        if (GrandCP(i,j)>0)
            ma(end+1) = GrandCP(i,j);
        end
    end
    AvgGrandCP(j) = mean(ma);
end

%% average only abs(TI)>0.05

for j = 1:size(GrandCP, 2)
    ma = [];
    for i = 1:size(GrandCP, 1)
        if(abs(TI(i))>0.05) & (GrandCP(i,j)>0)
            ma(end+1) = GrandCP(i,j);
        end
    end
    AvgGrandCPTI005(j) = mean(ma);
end


%%
figure(4565), clf, hold on,
plot([0:stepSize:20000]/10, AvgGrandCP);
plot([0:stepSize:20000]/10, AvgGrandCPTI005, 'r');
refline(0,0.5);


