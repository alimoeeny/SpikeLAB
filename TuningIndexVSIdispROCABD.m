% This file has been cross check with B. and is absolutely bug free as
% of 10/16/2009 for the AllABDNeurons.txt for the ABD and cylinder
% obviously.

clear;
clc

cd /b/ali/matlab/SpikeLAB/dev

try
    cd dev
catch tempE
    % we must be already inside dev
end
DataPath = GetDataPath('server');

ShowSingleCellSDFs = 0; % 0 or 1
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('TWO');
%Prep

% StartTime = 500;
% FinishTime = 10000;

doitsquare = 0;

% ABD
% load /Users/ali/DropBox/Projects/BCode/AllABDNeurons.mat
% AllNeurons = AllABDNeurons;
%  % AllNeurons = AllABDNeurons(57:end); disp ( ' O N L Y   I C A R U S ! ! ! !   B E W A R E ! ! !');
% clear AllABDNeurons;
% FileType = 'ABD';
% StimulusType = 'cylinder';
% StartTime  = 10000;  %10000; % 6500;
% FinishTime = 20000;

% % % DID
% load ../AllDIDNeurons.mat
% AllNeurons = AllDIDNeurons;
% clear AllDIDNeurons
% FileType = 'DID';
% StimulusType = 'cylinder';
%StartTime  = 10000; %10000; % 6500;
%FinishTime = 20000;


% % % DIDB
% load ../AllDIDBNeurons.mat
% AllNeurons = AllDIDBNeurons;
% clear AllDIDBNeurons
% FileType = 'DIDB';
% StimulusType = 'cylinder';
% StartTime  = 10000; %10000; % 6500;
% FinishTime = 20000;


% % TWO
% load('../AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWONeurons;
% FileType = 'TWO';
% StimulusType = 'cylinder';
% StartTime  = 500; %10000; %5500;
% FinishTime = 20000;

% % % DPI Cylinder
% load('../AllPursuitNeurons.mat');
% AllNeurons = AllPursuitNeurons;
% clear AllPursuitNeurons;
% FileType = 'DPI';
% StimulusType = 'cylinder';
% StartTime  = 500;% 10000; %500; %10000;
% FinishTime = 20000;


% % DPI rds
% load('../AllPursuitNeuronsrds.mat');
% AllNeurons = AllPursuitNeuronsrds;
% clear AllPursuitNeuronsrds;
% FileType = 'DPI';
% StimulusType = 'rds';
% StartTime  = 500;% 10000; %500; %10000;
% FinishTime = 20000;
%

% % JPI rds
% load('../AllPursuitNeuronsJPIrds.mat');
% AllNeurons = AllPursuitNeuronsJPIrds;
% clear AllPursuitNeuronsJPIrds;
% FileType = 'JPI';
% StimulusType = 'rds';
% StartTime  = 500;% 10000; %500; %10000;
% FinishTime = 20000;


% % % BDID
% load('../AllBDIDNeuronsALL.mat');
% AllNeurons = AllBDIDNeuronsALL;
% clear AllBDIDNeurons;
% FileType = 'BDID';
% StimulusType = 'cylinder';
% StartTime  = 10000;
% FinishTime = 20000;


%  AllNeurons =  SelectByMonkey(AllNeurons, 'ic');
%  AllNeurons =  SelectByMonkey(AllNeurons, 'dae');
%  disp('  O N E   M O N K E Y   A T  A  T I M E ');

%par
for iN=[1:length(AllNeurons)] %[1:33 40:length(AllNeurons)],
    if (sum(strcmpi(AllNeurons(iN), {'dae156','dae449', 'dae4492','dae5232', 'ic313'})>0))
        debug = 1;
    end
    
    IdColor{iN} = [1 0.5 0.1];
    DotSizes(iN) = 100;
    toolowFR = 0;
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber);
    TI(iN)      = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, [], DataPath);
    TIcyldx(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, 'DT',     [], DataPath);
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
    
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);
    %    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(filepath);
    Expt = Neuron.Expt;
    fileNames{iN} = filename;
    
    orS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 'or');
    backGroundOrS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 'bo');
        
    if strcmpi(StimulusType, 'cylinder')
        if(TI(iN)>0)
            pD = 1;
        elseif TI(iN)<0
            pD = 2;
        else
            pD = 0;
        end
        %pD      = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
        
        pDs(iN) = pD;
        
        if (pD == -1)
            disp('Ohoooy! pD? what are you D O I N G ! ? ');
        end
    else
        pD = -2;
    end
    
    if strcmpi(FileType, 'rds')
        prefRDSdx(iN) =  PreferredRDSDisparity(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    end
    %conditions = GetConditions(Expt, FileType, pD);
    try
        rdsPrefDirFailed = 0;
        rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName, DataPath);
    catch e
        rdsPrefDirFailed = 1;
        disp(['rdsPrefDir FAILED HERE ! ! ' , e.message]);
    end
    if (rdsPrefDirFailed || (length(rdsPrefDir)>1))
        rdsPrefDir = Expt.Stimvals.or;
    end
    prdsDir(iN) = rdsPrefDir;
    
    conditions = logical([]);
    
    if strcmpi(FileType, 'JPI')
        if isfield(Expt.Trials,'dfx')
            deltafxy = [Expt.Trials(:).dfx] - [Expt.Trials(:).fx];
            dfxy = [Expt.Trials(:).dfx];
        else
            deltafxy = [Expt.Trials(:).dfy] - [Expt.Trials(:).fy];
            dfxy = [Expt.Trials(:).dfy];
        end
        tempdelta = fix(deltafxy);
        tempd = fix(dfxy);
        for td = 1: length(tempdelta)
            if tempdelta(td) == 0
                tempdelta(td) = 0.5 * sign(deltafxy(td));
            end
        end
        deltafxy = tempdelta;
        Speeds = unique(abs(deltafxy));
        
        conditions(1,:) = deltafxy>0;
        conditions(2,:) = deltafxy<0;
        
        conditions(21,:) = deltafxy>0 & [Expt.Trials(:).jv]==0;
        conditions(22,:) = deltafxy<0 & [Expt.Trials(:).jv]==0;
        conditions(23,:) = deltafxy>0 & [Expt.Trials(:).jv]>0;
        conditions(24,:) = deltafxy<0 & [Expt.Trials(:).jv]>0;
        conditions(25,:) = deltafxy>0 & [Expt.Trials(:).jv]<0;
        conditions(26,:) = deltafxy<0 & [Expt.Trials(:).jv]<0;
        
        
    else if strcmpi(FileType, 'DPI')
            if strcmpi(StimulusType, 'cylinder')
                [conditions, r2p, r2n, Speeds, deltafxy] = GetConditions(Expt, FileType, pD, prdsDir(iN));
            else % rds
                [deltafxy, Speeds] = dpiDeltaFXY(Expt);
                conditions(1,:) = deltafxy>0;
                conditions(2,:) = deltafxy<0;
                if (isfield(Expt.Trials(1), 'dfx'))
                    pi = [Expt.Trials(:).dfx];
                else
                    pi = [Expt.Trials(:).dfy];
                end
                
                if ((Expt.Stimvals.or<=180) && (Expt.Stimvals.or>=-180))
                    conditions(20, :) = sign(pi) == sign(Expt.Stimvals.or);
                    conditions(21, :) = sign(pi) ~= sign(Expt.Stimvals.or);
                    
                    conditions(22, :) = (sign(pi) == sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) == sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0);
                    conditions(23, :) = (sign(pi) ~= sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) == sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0);
                    conditions(24, :) = (sign(pi) == sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) ~= sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0);
                    conditions(25, :) = (sign(pi) ~= sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) ~= sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0);
                    conditions(26, :) = sign(pi) == sign(Expt.Stimvals.or) & ([Expt.Trials(:).dx]==0);
                    conditions(27, :) = sign(pi) ~= sign(Expt.Stimvals.or) & ([Expt.Trials(:).dx]==0);
                else
                    conditions(20, :) = sign(pi) ~= sign(Expt.Stimvals.or);
                    conditions(21, :) = sign(pi) == sign(Expt.Stimvals.or);
                    
                    conditions(22, :) = ((sign(pi) ~= sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) == sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0));
                    conditions(23, :) = ((sign(pi) == sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) == sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0));
                    conditions(24, :) = ((sign(pi) ~= sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) ~= sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0));
                    conditions(25, :) = ((sign(pi) == sign(Expt.Stimvals.or)) & (([Expt.Trials(:).dx]>0) ~= sign(2 - prefRDSdx(iN))) & ([Expt.Trials(:).dx]~=0));
                    conditions(26, :) = sign(pi) ~= sign(Expt.Stimvals.or) & ([Expt.Trials(:).dx]==0);
                    conditions(27, :) = sign(pi) == sign(Expt.Stimvals.or) & ([Expt.Trials(:).dx]==0);
                end
            end
        else
            if strcmpi(FileType, 'BDID')
                if(mean([Expt.Trials([Expt.Trials(:).Id]>0).RespDir])>0)
                    ResponseToPositive = 1;
                    ResponseToNegative = -1;
                else
                    ResponseToPositive = -1;
                    ResponseToNegative = 1;
                end
            else
                if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
                    ResponseToPositive = 1;
                    ResponseToNegative = -1;
                else
                    ResponseToPositive = -1;
                    ResponseToNegative = 1;
                end
            end
            if strcmp(FileType, 'TWO')
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
                    %FLIP GRAND CP
                    conditions(12,:) = [Expt.Trials(:).RespDir]==ResponseToNegative & [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(13,:) = [Expt.Trials(:).RespDir]==ResponseToPositive & [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(14,:) = [Expt.Trials(:).RespDir]==ResponseToNegative & [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(15,:) = [Expt.Trials(:).RespDir]==ResponseToPositive & [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    
                    % 16 - 19 are the the two dxs near zero disparity.
                    conditions(16,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    conditions(17,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    conditions(18,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    conditions(19,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    
                    %these are very suspicious, why ~= ? what I was thinking
                    conditions(20,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToNegative;
                    conditions(21,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToPositive;
                    conditions(22,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToPositive;
                    conditions(23,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToNegative;
                    
                    
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
                    
                    %FLIP GRAND CP
                    conditions(12,:) = [Expt.Trials(:).RespDir]==ResponseToPositive & [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(13,:) = [Expt.Trials(:).RespDir]==ResponseToNegative & [Expt.Trials(:).bd]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(14,:) = [Expt.Trials(:).RespDir]==ResponseToPositive & [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    conditions(15,:) = [Expt.Trials(:).RespDir]==ResponseToNegative & [Expt.Trials(:).bd]== max([Expt.Trials(:).bd]) & [Expt.Trials(:).dx]== min([Expt.Trials(:).bd]) & [Expt.Trials(:).RespDir]~=0; 
                    
                    % 16 - 19 are the the two dxs near zero disparity.
                    conditions(16,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    conditions(17,:)= [Expt.Trials(:).dx] == min([Expt.Trials([Expt.Trials(:).dx]>0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    conditions(18,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    conditions(19,:)= [Expt.Trials(:).dx] == max([Expt.Trials([Expt.Trials(:).dx]<0).dx]) & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    
                    %these are very suspicious, why ~= ? what I was thinking
                    conditions(20,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToPositive;
                    conditions(21,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToNegative;
                    conditions(22,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToNegative;
                    conditions(23,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=ResponseToPositive;

                    % 24 - 27 are the the two dx zero disparitis for Z-Scored CP calculation.
                    conditions(24,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    conditions(25,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    conditions(26,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
                    conditions(27,:)= [Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
                    
                end
            else if strcmp(FileType, 'BDID')
                    conditions = GetConditions(Expt, FileType, pD);
                else
                    
                    % DID
                    conditions = GetConditions(Expt, FileType, pD);
                    
                end
            end
        end
    end
    if(sum(conditions(1,:))<7 || sum(conditions(2,:))<7)
        disp(strcat(num2str(iN) , ' - ', num2str(NeuronNumber, '%-4.3d') , ' - ', num2str(size([Expt.Trials(:)],1)), ' - TOO LOW TRIAL COUNT - BEWARE ! ! ! ', num2str(sum(conditions(1,:))), ' - ',num2str(sum(conditions(2,:)))));
        toolowFR = 1;
        continue;
    end
    
    
    if (FinishTime == 0)
        if(isfield(Expt.Trials, 'dur'))
            FinishTime = round(median([Expt.Trials(:).dur])) + 500;
        else
            FinishTime = round(median([Expt.Trials(:).End] - [Expt.Trials(:).Start])) + 500;
        end
    end
    SpikeCounts = zeros(length([Expt.Trials]),1);
    for tr = 1: length([Expt.Trials]),
        if doitsquare
            SpikeCounts(tr) = sqrt(sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime));
        else
            SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
        end
    end
    
    if (strfind(['dae494','ic343','ic377','ic107','ic311','ic367','','','','','','',], strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'))))
        debug = 1;
    end
    
    if (strcmpi(FileType, 'JPI'))
        dpJPIRDS1(iN) = dPrime(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
        dpJPIRDS2(iN) = dPrime(SpikeCounts(conditions(21,:)), SpikeCounts(conditions(22,:)));
        dpJPIRDS3(iN) = dPrime(SpikeCounts(conditions(23,:)), SpikeCounts(conditions(24,:)));
        dpJPIRDS4(iN) = dPrime(SpikeCounts(conditions(25,:)), SpikeCounts(conditions(26,:)));
    else if (strcmpi(FileType, 'DPI'))
            if strcmpi(StimulusType, 'cylinder')
                
                dprs = []; dprsdxp=[]; dprsdxn=[];
                for ss = 1: length(Speeds),
                    cs = (abs(deltafxy) == Speeds(ss));
                    dprs(ss) = dPrime(SpikeCounts(conditions(1,:) & cs), SpikeCounts(conditions(2,:) & cs));
                    dprsdxp(ss) = dPrime(SpikeCounts(conditions(3,:) & cs), SpikeCounts(conditions(4,:) & cs));
                    dprsdxn(ss) = dPrime(SpikeCounts(conditions(5,:) & cs), SpikeCounts(conditions(6,:) & cs)); % changed from 6 vs 5 to 5 vs 6 by Ali 8/18/11
                end
                dprs(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
                dprsdxp(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(4,:)));
                dprsdxn(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(5,:)), SpikeCounts(conditions(6,:))); % changed from 6 vs 5 to 5 vs 6 by Ali 8/18/11
                
                dprimes{iN} = [dprs dprsdxp dprsdxn];
                dpdxless(iN)= dPrime(SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:)));
                
                dpdx(iN)= dPrime(SpikeCounts(conditions(9,:)), SpikeCounts(conditions(10,:)));
                
                %        dppursdir(iN)= dPrime(SpikeCounts(conditions(11,:)), SpikeCounts(conditions(12,:)));
                
                
                PIFast(iN) = (mean(SpikeCounts(conditions(11,:))) - mean(SpikeCounts(conditions(12,:)))) ./ (mean(SpikeCounts(conditions(11,:))) + mean(SpikeCounts(conditions(12,:))));
                PISlow(iN) = (mean(SpikeCounts(conditions(21,:))) - mean(SpikeCounts(conditions(22,:)))) ./ (mean(SpikeCounts(conditions(21,:))) + mean(SpikeCounts(conditions(22,:))));
                PIMid(iN)  = (mean(SpikeCounts(conditions(31,:))) - mean(SpikeCounts(conditions(32,:)))) ./ (mean(SpikeCounts(conditions(31,:))) + mean(SpikeCounts(conditions(32,:))));
                
                PIOld(iN)  = (mean(SpikeCounts(conditions(7,:))) - mean(SpikeCounts(conditions(8,:)))) ./ (mean(SpikeCounts(conditions(7,:))) + mean(SpikeCounts(conditions(8,:))));
                PIOldZ(iN) = (mean(SpikeCounts(conditions(1,:))) - mean(SpikeCounts(conditions(2,:)))) ./ (mean(SpikeCounts(conditions(1,:))) + mean(SpikeCounts(conditions(2,:))));
                PIOldP(iN) = (mean(SpikeCounts(conditions(3,:))) - mean(SpikeCounts(conditions(4,:)))) ./ (mean(SpikeCounts(conditions(3,:))) + mean(SpikeCounts(conditions(4,:))));
                PIOldN(iN) = (mean(SpikeCounts(conditions(5,:))) - mean(SpikeCounts(conditions(6,:)))) ./ (mean(SpikeCounts(conditions(5,:))) + mean(SpikeCounts(conditions(6,:))));
                PI(iN) = (mean(SpikeCounts(conditions(1,:) | conditions(3,:) | conditions(5,:))) ...
                    - mean(SpikeCounts(conditions(2,:) | conditions(4,:) | conditions(6,:)))) ./ ...
                    (mean(SpikeCounts(conditions(1,:) | conditions(3,:) | conditions(5,:))) ...
                    + mean(SpikeCounts(conditions(2,:) | conditions(4,:) | conditions(6,:))));
                %         PIZero(iN)= (mean(SpikeCounts(conditions(11,:))) - mean(SpikeCounts(conditions(12,:)))) ./ (mean(SpikeCounts(conditions(11,:))) + mean(SpikeCounts(conditions(12,:))));
                %         PIPref(iN)= (mean(SpikeCounts(conditions(13,:))) - mean(SpikeCounts(conditions(14,:)))) ./ (mean(SpikeCounts(conditions(13,:))) + mean(SpikeCounts(conditions(14,:))));
                %         PINull(iN)= (mean(SpikeCounts(conditions(15,:))) - mean(SpikeCounts(conditions(16,:)))) ./ (mean(SpikeCounts(conditions(15,:))) + mean(SpikeCounts(conditions(16,:))));
                %         PI(iN) = mean([PIZero(iN), PIPref(iN), PINull(iN)]);
                %
                %         PIZerodx(iN) = (mean(SpikeCounts(conditions(22,:))) - mean(SpikeCounts(conditions(23,:)))) ./ (mean(SpikeCounts(conditions(22,:))) + mean(SpikeCounts(conditions(23,:))));
                %         PIPrefdx(iN) = (mean(SpikeCounts(conditions(24,:))) - mean(SpikeCounts(conditions(25,:)))) ./ (mean(SpikeCounts(conditions(24,:))) + mean(SpikeCounts(conditions(25,:))));
                %         PINulldx(iN) = (mean(SpikeCounts(conditions(26,:))) - mean(SpikeCounts(conditions(27,:)))) ./ (mean(SpikeCounts(conditions(26,:))) + mean(SpikeCounts(conditions(27,:))));
                %         PINew(iN) = mean([PIZerodx(iN), PIPrefdx(iN), PINulldx(iN)]);
                
                NS{iN} = {SpikeCounts, conditions, deltafxy};
            else
                dpRDS1(iN) = dPrime(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
                dpRDS2(iN) = dPrime(SpikeCounts(conditions(21,:)), SpikeCounts(conditions(22,:)));
                dpRDS3(iN) = dPrime(SpikeCounts(conditions(23,:)), SpikeCounts(conditions(24,:)));
                dpRDS4(iN) = dPrime(SpikeCounts(conditions(25,:)), SpikeCounts(conditions(26,:)));
                
                %           PIOld(iN) = (mean(SpikeCounts(conditions(7,:))) - mean(SpikeCounts(conditions(8,:)))) ./ (mean(SpikeCounts(conditions(7,:))) + mean(SpikeCounts(conditions(8,:))));
                %         PIOldZ(iN) = (mean(SpikeCounts(conditions(1,:))) - mean(SpikeCounts(conditions(2,:)))) ./ (mean(SpikeCounts(conditions(1,:))) + mean(SpikeCounts(conditions(2,:))));
                %         PIOldP(iN) = (mean(SpikeCounts(conditions(3,:))) - mean(SpikeCounts(conditions(4,:)))) ./ (mean(SpikeCounts(conditions(3,:))) + mean(SpikeCounts(conditions(4,:))));
                %         PIOldN(iN) = (mean(SpikeCounts(conditions(5,:))) - mean(SpikeCounts(conditions(6,:)))) ./ (mean(SpikeCounts(conditions(5,:))) + mean(SpikeCounts(conditions(6,:))));
                %
                %         PIZero(iN)= (mean(SpikeCounts(conditions(1,:))) - mean(SpikeCounts(conditions(2,:)))) ./ (mean(SpikeCounts(conditions(1,:))) + mean(SpikeCounts(conditions(2,:))));
                %           PIPref(iN)= (mean(SpikeCounts(conditions(3,:))) - mean(SpikeCounts(conditions(4,:)))) ./ (mean(SpikeCounts(conditions(3,:))) + mean(SpikeCounts(conditions(4,:))));
                %           PINull(iN)= (mean(SpikeCounts(conditions(5,:))) - mean(SpikeCounts(conditions(6,:)))) ./ (mean(SpikeCounts(conditions(5,:))) + mean(SpikeCounts(conditions(6,:))));
                %           PI(iN) = mean([PIZero(iN), PIPref(iN), PINull(iN)]);
                %
                %           PIZerodx(iN) = (mean(SpikeCounts(conditions(22,:))) - mean(SpikeCounts(conditions(23,:)))) ./ (mean(SpikeCounts(conditions(22,:))) + mean(SpikeCounts(conditions(23,:))));
                %           PIPrefdx(iN) = (mean(SpikeCounts(conditions(24,:))) - mean(SpikeCounts(conditions(25,:)))) ./ (mean(SpikeCounts(conditions(24,:))) + mean(SpikeCounts(conditions(25,:))));
                %           PINulldx(iN) = (mean(SpikeCounts(conditions(26,:))) - mean(SpikeCounts(conditions(27,:)))) ./ (mean(SpikeCounts(conditions(26,:))) + mean(SpikeCounts(conditions(27,:))));
                %           PINew(iN) = mean([PIZerodx(iN), PIPrefdx(iN), PINulldx(iN)]);
            end
        else
            
            if (strcmpi(FileType, 'BDID'))
                %IdBiasROC1(iN) = ROCAUC(SpikeCounts(conditions(3,:)|conditions(4,:)), SpikeCounts(conditions(5,:)|conditions(6,:)));
                IdBiasROC1(iN) = ROCAUC(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(5,:)|conditions(6,:)));
                [rr, pp] = ROCAUCSignificance(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(5,:)|conditions(6,:)));
                IdBiasROCSig(iN) = pp;
                
                IdBiasROC2(iN) = ROCAUC(SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:)));
                [rr, pp] = ROCAUCSignificance(SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:)));
                IdBiasROCSig2(iN) = pp;
                ROCpairs1{iN} = {SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:))};
                
                
                IdBiasROC3(iN) = ROCAUC(SpikeCounts(conditions(9,:)), SpikeCounts(conditions(10,:)));
                [rr, pp] = ROCAUCSignificance(SpikeCounts(conditions(9,:)), SpikeCounts(conditions(10,:)));
                IdBiasROCSig3(iN) = pp;
                
                %         bdCrossTalk(iN) = 0.5 * abs( (sum(conditions(16,:)) - sum(conditions(17,:))) / (sum(conditions(16,:)) + sum(conditions(17,:))) + ...
                %                                      (sum(conditions(18,:)) - sum(conditions(19,:))) / (sum(conditions(18,:)) + sum(conditions(19,:))) );
                %
                %         bdCrossTalk2(iN)= 0.5 * abs( (sum(conditions(16,:)) - sum(conditions(19,:))) / (sum(conditions(16,:)) + sum(conditions(19,:))) + ...
                %                                      (sum(conditions(18,:)) - sum(conditions(17,:))) / (sum(conditions(18,:)) + sum(conditions(17,:))) );
                %
                %         bdCrossTalk3(iN)= 0.5 *    ( (sum(conditions(16,:)) - sum(conditions(19,:))) / (sum(conditions(16,:)) + sum(conditions(19,:))) + ...
                %                                      (sum(conditions(18,:)) - sum(conditions(17,:))) / (sum(conditions(18,:)) + sum(conditions(17,:))) );
                %
                %         bdCrossTalk4(iN)= 0.5 * abs( (sum(conditions(16,:)) - sum(conditions(19,:))) / (sum(conditions(16,:)) + sum(conditions(19,:))))+ ...
                %                                 abs( (sum(conditions(18,:)) - sum(conditions(17,:))) / (sum(conditions(18,:)) + sum(conditions(17,:))));
                %
                
                boS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 'bo');
                if (orS(iN)==90)
                    c1 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                    c2 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                    c3 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                    c4 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                elseif (orS(iN)==-90)
                    c1 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                    c2 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                    c3 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                    c4 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                elseif (orS(iN)==0)
                    c1 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                    c2 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                    c3 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                    c4 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                elseif (orS(iN)==180)
                    c1 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                    c2 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                    c3 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                    c4 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                else
                    c1 = 0; c2 = 0; c3 = 0; c4 = 0;
                    disp(['ORIENTATION IS :  ', num2str(orS(iN)), '  WWHHAATT CCAANN II DDOO == == == == == == == =='])
                end
                bdCrossTalk5(iN) = 0.5 * ( ...
                    ((sum(c2) - sum(c1)) / (sum(c1) + sum(c2))) + ...
                    ((sum(c3) - sum(c4)) / (sum(c4) + sum(c3))) );
                
                %          bdCrossTalk6(iN) = 0.5 * ( ...
                %                             ((sum(c1) - sum(c2)) / (sum(c1) + sum(c2))) + ...
                %                             ((sum(c3) - sum(c4)) / (sum(c4) + sum(c3))) );
                
                % Neuronal effect of Id for large bd vs small bd trails
                lbdt = abs([Expt.Trials(:).bd]) == max(abs([Expt.Trials(:).bd]));
                IdBiasROC2LargebdTrials(iN) = ROCAUC(SpikeCounts(conditions(7,:) & lbdt), SpikeCounts(conditions(8,:) & lbdt));
                sbdt = abs([Expt.Trials(:).bd]) < 0.01;
                IdBiasROC2SmallbdTrials(iN) = ROCAUC(SpikeCounts(conditions(7,:) & sbdt), SpikeCounts(conditions(8,:) & sbdt));
                
                % Choice Probability for Id trials
                IdCP1(iN) = ROCAUC(SpikeCounts(conditions(16,:)), SpikeCounts(conditions(17,:)));
                IdCP2(iN) = ROCAUC(SpikeCounts(conditions(19,:)), SpikeCounts(conditions(18,:)));
                IdCP3(iN) = 0.5 * (IdCP1(iN) + IdCP2(iN));
                IdCP4(iN) = ROCAUC(SpikeCounts(conditions(16,:) | conditions(19,:)), SpikeCounts(conditions(17,:) | conditions(18,:)));
                
            else % DID, ABD, TWO, DIDB, ...
                
                mtn = 10; % minimum acceptable trials count
                if ((sum(conditions(1,:))>=mtn) && ((sum(conditions(2,:))>=mtn) )) ,
                    IdBiasROC1(iN) = ROCAUC(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
                    [rr, pp] = ROCAUCSignificance(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
                    IdBiasROCSig(iN) = pp;
                else
                    IdBiasROC1(iN) = -10;
                    IdBiasROCSig(iN) = -10;
                end
                ROCpairs1{iN} = {SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:))};
                
                if(~strcmpi(FileType, 'TWO'))
                    Next2ZeroROC1(iN) = ROCAUC(SpikeCounts(conditions(13,:)), SpikeCounts(conditions(14,:)));
                    % next to zero z scored
                    
                    a = zscore([SpikeCounts(conditions(16,:)); SpikeCounts(conditions(17,:))]);
                    b = zscore([SpikeCounts(conditions(18,:)); SpikeCounts(conditions(19,:))]);
                    aa = [a(1 : sum(conditions(16,:))) ; b(1 : sum(conditions(18,:)))];
                    bb = [a(sum(conditions(16,:))+1 : sum(conditions(16,:))+sum(conditions(17,:))) ; b(sum(conditions(18,:))+1 : sum(conditions(18,:))+sum(conditions(19,:)))];
                    Next2ZeroROCZScored(iN) = ROCAUC(aa, bb);
                    Next2ZeroROCPref(iN) = ROCAUC(SpikeCounts(conditions(16,:)), SpikeCounts(conditions(17,:)));
                    Next2ZeroROCNull(iN) = ROCAUC(SpikeCounts(conditions(18,:)), SpikeCounts(conditions(19,:)));
                    w1 = min(sum(conditions(16,:)),sum(conditions(17,:))); if (w1==0), w1 = max(sum(conditions(16,:)),sum(conditions(17,:))); end
                    w2 = min(sum(conditions(18,:)),sum(conditions(19,:))); if (w2==0), w2 = max(sum(conditions(18,:)),sum(conditions(19,:))); end
                    if (Next2ZeroROCPref(iN) == -1)
                        Next2ZeroROCPNWeigh(iN) = Next2ZeroROCNull(iN);
                    else
                        Next2ZeroROCPNWeigh(iN) = Next2ZeroROCPref(iN);
                    end
                    if ((Next2ZeroROCPref(iN) ~= -1) && (Next2ZeroROCNull(iN) ~= -1))
                        Next2ZeroROCPNWeigh(iN)= ((Next2ZeroROCPref(iN) * w1) + (Next2ZeroROCNull(iN) * w2)) / (w1 + w2);
                    end
                    
                    
                    
                    % Z-scored CP
                    a = zscore([SpikeCounts(conditions(24,:)); SpikeCounts(conditions(25,:))]);
                    b = zscore([SpikeCounts(conditions(26,:)); SpikeCounts(conditions(27,:))]);
                    aa = [a(1 : sum(conditions(24,:))) ; b(1 : sum(conditions(26,:)))];
                    bb = [a(sum(conditions(24,:))+1 : sum(conditions(24,:))+sum(conditions(25,:))) ; b(sum(conditions(26,:))+1 : sum(conditions(26,:))+sum(conditions(27,:)))];
                    CPatZeroROCZScored(iN) = ROCAUC(aa, bb);
                    CPatZeroROCPref(iN) = ROCAUC(SpikeCounts(conditions(24,:)), SpikeCounts(conditions(25,:)));
                    CPatZeroROCNull(iN) = ROCAUC(SpikeCounts(conditions(26,:)), SpikeCounts(conditions(27,:)));
                    w1 = min(sum(conditions(24,:)),sum(conditions(25,:))); if (w1==0), w1 = max(sum(conditions(24,:)),sum(conditions(25,:))); end
                    w2 = min(sum(conditions(26,:)),sum(conditions(27,:))); if (w2==0), w2 = max(sum(conditions(26,:)),sum(conditions(27,:))); end
                    if (CPatZeroROCPref(iN) == -1)
                        CPatZeroROCPNWeigh(iN) = CPatZeroROCNull(iN);
                    else
                        CPatZeroROCPNWeigh(iN) = CPatZeroROCPref(iN);
                    end
                    if ((CPatZeroROCPref(iN) ~= -1) && (CPatZeroROCNull(iN) ~= -1))
                        CPatZeroROCPNWeigh(iN)= ((CPatZeroROCPref(iN) * w1) + (CPatZeroROCNull(iN) * w2)) / (w1 + w2);
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
                    aaNZ = []; bbNZ = []; AANZ={}; BBNZ={}; %Grand CP winthout zero dx trials
                    %CP in next to zero dx trials with null DX (no id here)
                    if ((sum(conditions(16,:))>=mtn) && ((sum(conditions(17,:))>=mtn) )) ,
                        at = pdx(1:sum(conditions(16,:)));
                        bt = pdx(sum(conditions(16,:))+1:sum(conditions(16,:))+sum(conditions(17,:)));
                        aa = [aa; at];
                        bb = [bb; bt];
                        AA{end+1} = at;
                        BB{end+1} = bt;
                        aaNZ = [aaNZ; at];
                        bbNZ = [bbNZ; bt];
                        AANZ{end+1} = at;
                        BBNZ{end+1} = bt;
                    end
                    %CP in next to zero dx trials with pref DX (no id here)
                    if ((sum(conditions(18,:))>=mtn) && ((sum(conditions(19,:))>=mtn) )) ,
                        at = ndx(1:sum(conditions(18,:)));
                        bt = ndx(sum(conditions(18,:))+1:sum(conditions(18,:))+sum(conditions(19,:)));
                        aa = [aa; at];
                        bb = [bb; bt];
                        AA{end+1} = at;
                        BB{end+1} = bt;
                        aaNZ = [aaNZ; at];
                        bbNZ = [bbNZ; bt];
                        AANZ{end+1} = at;
                        BBNZ{end+1} = bt;
                    end
                    
                    GrandCP(iN) = ROCAUC(aa, bb);
                    [roc2compare, p, pp] = ROCAUCSignificanceGrandCP(AA, BB);
                    GrandCPSig(iN) = p;
                    GrandCPNZ(iN) = ROCAUC(aaNZ, bbNZ);
                    
                    if(GrandCP(iN)<0)
                        debug = 1;
                    end
                else % TWO <- GrandCP and more ----------------------------------------------------------------------------------------------------------------------------------
                    
                    % next to zero z scored
                    a = zscore([SpikeCounts(conditions(16,:)); SpikeCounts(conditions(17,:))]);
                    b = zscore([SpikeCounts(conditions(18,:)); SpikeCounts(conditions(19,:))]);
                    aa = [a(1 : sum(conditions(16,:))) ; b(1 : sum(conditions(18,:)))];
                    bb = [a(sum(conditions(16,:))+1 : sum(conditions(16,:))+sum(conditions(17,:))) ; b(sum(conditions(18,:))+1 : sum(conditions(18,:))+sum(conditions(19,:)))];
                    Next2ZeroROCZScored(iN) = ROCAUC(aa, bb);
                    Next2ZeroROCPref(iN) = ROCAUC(SpikeCounts(conditions(16,:)), SpikeCounts(conditions(17,:)));
                    Next2ZeroROCNull(iN) = ROCAUC(SpikeCounts(conditions(18,:)), SpikeCounts(conditions(19,:)));
                    w1 = min(sum(conditions(16,:)),sum(conditions(17,:))); if (w1==0), w1 = max(sum(conditions(16,:)),sum(conditions(17,:))); end
                    w2 = min(sum(conditions(18,:)),sum(conditions(19,:))); if (w2==0), w2 = max(sum(conditions(18,:)),sum(conditions(19,:))); end
                    if (Next2ZeroROCPref(iN) == -1)
                        Next2ZeroROCPNWeigh(iN) = Next2ZeroROCNull(iN);
                    else
                        Next2ZeroROCPNWeigh(iN) = Next2ZeroROCPref(iN);
                    end
                    if ((Next2ZeroROCPref(iN) ~= -1) && (Next2ZeroROCNull(iN) ~= -1))
                        Next2ZeroROCPNWeigh(iN)= ((Next2ZeroROCPref(iN) * w1) + (Next2ZeroROCNull(iN) * w2)) / (w1 + w2);
                    end
                    
                    
                    
                    % Z-scored CP
                    a = zscore([SpikeCounts(conditions(24,:)); SpikeCounts(conditions(25,:))]);
                    b = zscore([SpikeCounts(conditions(26,:)); SpikeCounts(conditions(27,:))]);
                    aa = [a(1 : sum(conditions(24,:))) ; b(1 : sum(conditions(26,:)))];
                    bb = [a(sum(conditions(24,:))+1 : sum(conditions(24,:))+sum(conditions(25,:))) ; b(sum(conditions(26,:))+1 : sum(conditions(26,:))+sum(conditions(27,:)))];
                    CPatZeroROCZScored(iN) = ROCAUC(aa, bb);
                    CPatZeroROCPref(iN) = ROCAUC(SpikeCounts(conditions(24,:)), SpikeCounts(conditions(25,:)));
                    CPatZeroROCNull(iN) = ROCAUC(SpikeCounts(conditions(26,:)), SpikeCounts(conditions(27,:)));
                    w1 = min(sum(conditions(24,:)),sum(conditions(25,:))); if (w1==0), w1 = max(sum(conditions(24,:)),sum(conditions(25,:))); end
                    w2 = min(sum(conditions(26,:)),sum(conditions(27,:))); if (w2==0), w2 = max(sum(conditions(26,:)),sum(conditions(27,:))); end
                    if (CPatZeroROCPref(iN) == -1)
                        CPatZeroROCPNWeigh(iN) = CPatZeroROCNull(iN);
                    else
                        CPatZeroROCPNWeigh(iN) = CPatZeroROCPref(iN);
                    end
                    if ((CPatZeroROCPref(iN) ~= -1) && (CPatZeroROCNull(iN) ~= -1))
                        CPatZeroROCPNWeigh(iN)= ((CPatZeroROCPref(iN) * w1) + (CPatZeroROCNull(iN) * w2)) / (w1 + w2);
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
                    aaNZ = []; bbNZ = []; AANZ={}; BBNZ={}; %Grand CP winthout zero dx trials
                    %CP in next to zero dx trials with null DX (no id here)
                    if ((sum(conditions(16,:))>=mtn) && ((sum(conditions(17,:))>=mtn) )) ,
                        at = pdx(1:sum(conditions(16,:)));
                        bt = pdx(sum(conditions(16,:))+1:sum(conditions(16,:))+sum(conditions(17,:)));
                        aa = [aa; at];
                        bb = [bb; bt];
                        AA{end+1} = at;
                        BB{end+1} = bt;
                        aaNZ = [aaNZ; at];
                        bbNZ = [bbNZ; bt];
                        AANZ{end+1} = at;
                        BBNZ{end+1} = bt;
                    end
                    %CP in next to zero dx trials with pref DX (no id here)
                    if ((sum(conditions(18,:))>=mtn) && ((sum(conditions(19,:))>=mtn) )) ,
                        at = ndx(1:sum(conditions(18,:)));
                        bt = ndx(sum(conditions(18,:))+1:sum(conditions(18,:))+sum(conditions(19,:)));
                        aa = [aa; at];
                        bb = [bb; bt];
                        AA{end+1} = at;
                        BB{end+1} = bt;
                        aaNZ = [aaNZ; at];
                        bbNZ = [bbNZ; bt];
                        AANZ{end+1} = at;
                        BBNZ{end+1} = bt;                    
                    end
                    
                    GrandCP(iN) = ROCAUC(aa, bb);
                    [roc2compare, p, pp] = ROCAUCSignificanceGrandCP(AA, BB);
                    GrandCPSig(iN) = p;
                    
                    GrandCPNZ(iN) = ROCAUC(aaNZ, bbNZ);
                    
                    if(GrandCP(iN)<0)
                        debug = 1;
                    end
                    
                    %% Contingency Tables
                    CT = zeros(2,2,2); % 3rd two, one for zscored rates and one for trial counts
                    CT(1,1,2) =  sum(conditions(24,:)); % pref bd pref choice
                    CT(1,2,2) =  sum(conditions(25,:)); % pref bd null choice
                    CT(2,2,2) =  sum(conditions(27,:)); % null bd null choice
                    CT(2,1,2) =  sum(conditions(26,:)); % null bd pref choice
                    
                    
                    zsc = zscore([SpikeCounts(conditions(24,:)); SpikeCounts(conditions(25,:));SpikeCounts(conditions(27,:));SpikeCounts(conditions(26,:))]);
                    bi = 1; ei = sum(conditions(24,:));
                    CT(1,1,1) =  mean(zsc(bi:ei));  % pref bd pref choice
                    bi = ei + 1; ei = ei + sum(conditions(25,:));
                    CT(1,2,1) =  mean(zsc(bi:ei));  % pref bd null choice
                    bi = ei + 1; ei = ei + sum(conditions(27,:));
                    CT(2,2,1) =  mean(zsc(bi:ei)); % null bd null choice
                    bi = ei + 1; ei = ei + sum(conditions(26,:));
                    CT(2,1,1) =  mean(zsc(bi:ei)); % null bd pref choice
                    
                    if (sum(isnan(CT(:)))>0)
                        debug = 1;
                    end
                    
                    ContignGTables{iN} = CT;
                    
                    
                    %% GRAND CP FOR FLIP
                    zpf = zscore([SpikeCounts(conditions(12,:)); SpikeCounts(conditions(13,:))]);
                    znf = zscore([SpikeCounts(conditions(14,:)); SpikeCounts(conditions(15,:))]);
                    aaF = []; bbF = []; AAF={}; BBF={};
                    % FLIP trials with null id
                    if ((sum(conditions(12,:))>=mtn) && ((sum(conditions(13,:))>=mtn) )) ,
                        at = zpf(1:sum(conditions(12,:)));
                        bt = zpf(sum(conditions(12,:))+1:sum(conditions(12,:))+sum(conditions(13,:)));
                        aaF = [aaF; at];
                        bbF = [bbF; bt];
                        AAF{end+1} = at;
                        BBF{end+1} = bt;
                    end
                    % FLIP trials with pref id
                    if ((sum(conditions(14,:))>=mtn) && ((sum(conditions(15,:))>=mtn) )) ,
                        at = znf(1:sum(conditions(14,:)));
                        bt = znf(sum(conditions(14,:))+1:sum(conditions(14,:))+sum(conditions(15,:)));
                        aaF = [aaF; at];
                        bbF = [bbF; bt];
                        AAF{end+1} = at;
                        BBF{end+1} = bt;
                    end
                    GrandCPFlip(iN) = ROCAUC(aaF, bbF);
                    
                    %% CP only with Biased trials
                    
                    %somthing is very suspicious about conditions(20 through 23 make sure that is what you mean here
                    CPOnlyCorrectBiased(iN) = ROCAUC(SpikeCounts(conditions(20,:)), SpikeCounts(conditions(22,:)));
                    CPOnlyWrongBiased(iN) = ROCAUC(SpikeCounts(conditions(21,:)), SpikeCounts(conditions(23,:)));
                end
                
                orS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
                
            end
            switch FileType
                case {'TWO', 'BDID'}
                    if Expt.Stimvals.bo == Expt.Stimvals.or
                        Manip(iN) = 1;
                    else
                        Manip(iN) = 2;
                    end
                otherwise
                    %CP(iN)      = ROCAUC(SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:)));
                    CP(iN)      = ROCAUC(SpikeCounts(conditions(20,:)), SpikeCounts(conditions(21,:)));
                    CPPref(iN)  = ROCAUC(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(5,:))); %CPPref(iN)  = ROCAUC(SpikeCounts(conditions(9,:)), SpikeCounts(conditions(10,:)));
                    CPNull(iN)  = ROCAUC(SpikeCounts(conditions(4,:)), SpikeCounts(conditions(6,:))); %CPNull(iN)  = ROCAUC(SpikeCounts(conditions(11,:)), SpikeCounts(conditions(12,:)));
                    %CP1(iN)     = ROCAUC(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(5,:)));
                    %CP2(iN)     = ROCAUC(SpikeCounts(conditions(6,:)), SpikeCounts(conditions(4,:)));
                    %CPm(iN)     = mean([CP1(iN) CP2(iN)]);
                    %CPeffect(iN)= CP1(iN) - CP2(iN);
            end
            if strcmp(FileType,'TWO')
                if (sum([63 469 130]==NeuronNumber)>0)
                    debug = 1;
                end
                if((sum(conditions(7,:)) + sum(conditions(8,:)) == 0) || (sum(conditions(10,:)) + sum(conditions(11,:)) == 0))
                    debug = 1;
                    disp(['missed them #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ #$ ', AllNeurons{iN}]);
                end
                mtn = 0; % minimum trials needed
                fa = -1; fb = -1;
                if (sum(conditions(8,:)) > mtn && sum(conditions(7,:))>mtn)
                    fa = ROCAUC(SpikeCounts(conditions(8,:)), SpikeCounts(conditions(7,:)));
                end
                if (sum(conditions(10,:)) > mtn && sum(conditions(11,:))>mtn)
                    fb = ROCAUC(SpikeCounts(conditions(10,:)), SpikeCounts(conditions(11,:)));
                end
                
                if abs(fa * fb ) ~= (fa * fb)
                    disp('Fishy - - - - - - - - - - - 9 9 9 9 9 9 8  8 8 8( ( ((  ** * * * ');
                    fa = 0; fb = 0;
                end
                if (fa > 0 && fb > 0)
                    FlipROC(iN) = 0.5 * ( fa + fb );
                else
                    if (fa>0)
                        FlipROC(iN) = fa;
                    else
                        FlipROC(iN) = fb;
                    end
                end
                ROCpairsFlip{iN} = {SpikeCounts(conditions(8,:)), SpikeCounts(conditions(7,:))};
                ROCpairsFlip{iN+ length(AllNeurons)} = {SpikeCounts(conditions(10,:)), SpikeCounts(conditions(11,:))};
            end
            
            % Bias Effectiveness
            if ~strcmp(FileType, 'BDID')
                if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
                    ResponseToPositive = 1;
                    ResponseToNegative = -1;
                else
                    ResponseToPositive = -1;
                    ResponseToNegative = 1;
                end
                
                if strcmp(FileType,'TWO')
                    BiaEff(iN) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir] ~= 0) ...
                        - 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir] ~= 0);
                    BiaEffFlip(iN) = 0.5 * ( ...
                        100 * sum([Expt.Trials(:).dx] == max([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] == [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == max([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] == [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] ~= 0) ...
                        - 100 * sum([Expt.Trials(:).dx] == max([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] ~= [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == max([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] ~= [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] ~= 0) ...
                        + 100 * sum([Expt.Trials(:).dx] == min([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] == [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == min([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] == [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] ~= 0) ...
                        - 100 * sum([Expt.Trials(:).dx] == min([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] ~= [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == min([Expt.Trials(:).dx]) & [Expt.Trials(:).bd] ~= [Expt.Trials(:).dx] & [Expt.Trials(:).RespDir] ~= 0) );
                    if isnan(BiaEffFlip(iN))
                        debug = 1;
                    end
                    xo = Expt.Stimvals.xo;
                    yo = Expt.Stimvals.yo;
                    if (0<sum(NeuronNumber==[ 130 128 125 324 330]))
                        RFproximity(iN) = -1;
                        FxProximity(iN) = sqrt(xo^2 + yo^2); % - (szrf);
                        disp(['----------------------------------------', num2str([FxProximity(iN),xo, yo])]);
                    else
                        if (isfield(Expt.Trials,'backxo'))
                            bx = median([Expt.Trials(:).backxo]);
                        else
                            bx = median([Expt.Stimvals.backxo]);
                        end
                        if (isfield(Expt.Trials,'backyo'))
                            by = median([Expt.Trials(:).backyo]);
                        else
                            by = median([Expt.Stimvals.backyo]);
                        end
                        sz = median([Expt.Trials(:).sz]);
                        rf = RFFit(MonkeyName, NeuronNumber, ClusterName, 0);
                        szrf = rf(1) + sz;
                        RFproximity(iN) = sqrt((xo - bx )^2 + (yo - by) ^2 ) - (szrf);
                        FxProximity(iN) = sqrt(xo^2 + yo^2); % - (szrf);
                        disp(['----------------------------------------', num2str([FxProximity(iN), xo, yo])]);
                    end
                else
                    BiaEff(iN) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0) ...
                        - 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
                end
            end
        end
    end
    if (strcmp(FileType, 'BDID') & abs(Expt.Stimvals.bo - Expt.Stimvals.or)~=90)% & (  Expt.Stimvals.or == 0 | Expt.Stimvals.or == 180 | Expt.Stimvals.or == -180))
        disp('Are you kidding! this is a * TWO Experiment');
        disp([IdBiasROC1(iN) Expt.Stimvals.or Expt.Stimvals.bo ]);
        %         disp(Expt.Stimvals.xo);
        %         disp(Expt.Stimvals.yo);
        thisisTWO(iN) = 1;
        debug = 1;
    end
    
    if(TI(iN)>0.077)
        IdColor{iN} = [1 0 0];
        DotSizes(iN) = 70;
    else if (TI(iN)<-0.077)
            IdColor{iN} = [0 0 1];
            DotSizes(iN) = 70;
        else if(TI(iN) < 0)
                IdColor{iN} = [0.3 0.6 1];
                DotSizes(iN) = 40;
            else if(TI(iN) > 0)
                    IdColor{iN} = [1 0.6 0.3];
                    DotSizes(iN) = 40;
                else if (TI(iN) == 0 || isnan(TI(iN)))
                        IdColor{iN} = [0 0 0];
                        DotSizes(iN) = 100;
                    end
                end
            end
        end
    end
    if strcmpi(MonkeyName, 'icarus')
        IdColor{iN} = [0 0 0];
    end
    if toolowFR == 1
        DotSizes(iN) = 120;
    else
        DotSizes(iN) = 40;
    end
    
    binoc{iN} = Expt.Stimvals.ve;
end

for iN = 1:length(AllNeurons)
    if(AllNeurons{iN}(1) == 'i'),
        monkeyTag(iN) = 1;
    end,
end

pDsT = (((pDs/2-1)*-2)*2)-1;

%%
save(['~/Desktop/matlab.', [FileType], '.', num2str(now) '.mat'], '-v7.3')

shart1 = pDsT==sign(TI) ;
shart2 = abs(TI)>0.05;
shart = shart1 & shart2;

%%
if strcmpi(FileType, 'DID')
    figure(1112), clf, hold on,
    clickscatter(TI, CPatZeroROCPNWeigh,  8, 8, fileNames);
    clickscatter(TI, Next2ZeroROCPNWeigh, 7, 8, fileNames);
    xlabel('Tuning Index');
    ylabel('CP (Green is CP at zero disparity and \r\nRed is CP right next to zero disparity)');
    refline(0.0, 0.5);
    
    figure(1122), clf, hold on,
    clickscatter(TI, CPatZeroROCZScored,  8, 8, fileNames);
    clickscatter(TI, Next2ZeroROCZScored, 7, 8, fileNames);
    xlabel('Tuning Index');
    ylabel('CP (Green is CP at zero disparity and \r\nRed is CP right next to zero disparity)');
    refline(0.0, 0.5);
    
    figure(3332), clf, hold on,
    clickscatter(Next2ZeroROCPNWeigh, CPatZeroROCPNWeigh, 6, 8, fileNames);
    clickscatter(Next2ZeroROCPNWeigh(abs(TI)>0.1), CPatZeroROCPNWeigh(abs(TI)>0.1), 5, 8, fileNames(abs(TI)>0.1));
    clickscatter(Next2ZeroROCPNWeigh(pDs==2), CPatZeroROCPNWeigh(pDs==2), 2, 8, fileNames(pDs==2));
    clickscatter(Next2ZeroROCPNWeigh(pDs==1), CPatZeroROCPNWeigh(pDs==1), 1, 8, fileNames(pDs==1));
    xlabel('CP right next to zero disparity');
    ylabel('CP at zero disparity');
    refline(0.0, 0.5);
    reflinexy(0.5, 1);
    
    figure(3322), clf, hold on,
    clickscatter(Next2ZeroROCZScored, CPatZeroROCZScored, 6, 8, fileNames);
    clickscatter(Next2ZeroROCZScored(abs(TI)>0.1), CPatZeroROCZScored(abs(TI)>0.1), 5, 8, fileNames(abs(TI)>0.1));
    clickscatter(Next2ZeroROCZScored(pDs==2), CPatZeroROCZScored(pDs==2), 2, 8, fileNames(pDs==2));
    clickscatter(Next2ZeroROCZScored(pDs==1), CPatZeroROCZScored(pDs==1), 1, 8, fileNames(pDs==1));
    xlabel('CP right next to zero disparity');
    ylabel('CP at zero disparity');
    refline(0.0, 0.5);
    reflinexy(0.5, 1);
    
    figure(4442), clf, hold on,
    clickscatter(Next2ZeroROCPNWeigh, Next2ZeroROCZScored, 3, 8, fileNames);
    clickscatter(Next2ZeroROCPNWeigh(abs(TI)>0.1), Next2ZeroROCZScored(abs(TI)>0.1), 4, 8, fileNames(abs(TI)>0.1));
    xlabel('CP right next to zero disparity - Weighted avg');
    ylabel('CP right next to zero disparity - ZScored');
    refline(0.0, 0.5);
    reflinexy(0.5, 1);
    
    figure(4422), clf, hold on,
    clickscatter(CPatZeroROCPNWeigh, CPatZeroROCZScored, 3, 8, fileNames);
    clickscatter(CPatZeroROCPNWeigh(abs(TI)>0.1), CPatZeroROCZScored(abs(TI)>0.1), 4, 8, fileNames(abs(TI)>0.1));
    xlabel('CP at zero disparity - Weighted avg');
    ylabel('CP at zero disparity - ZScored');
    refline(0.0, 0.5);
    reflinexy(0.5, 1);
    %%
    figure(9999), clf, hold on
    %clickscatter(GrandCP, IdBiasROC1, 8, 8, fileNames);
    %clickscatter(GrandCP(abs(TI)>0.1), IdBiasROC1(abs(TI)>0.1), 6, 8, fileNames(abs(TI)>0.1));
    clickscatter(GrandCP(pDs==2), IdBiasROC1(pDs==2), 5, 8, fileNames(pDs==2));
    clickscatter(GrandCP(pDs==1), IdBiasROC1(pDs==1), 4, 8, fileNames(pDs==1));
    %clickscatter(GrandCP((pDs==2)&(abs(TI)>0.1)), IdBiasROC1((pDs==2)&(abs(TI)>0.1)), 4, 8, fileNames((pDs==2) & (abs(TI)>0.1)));
    %clickscatter(GrandCP((pDs==1)&(abs(TI)>0.1)), IdBiasROC1((pDs==1)&(abs(TI)>0.1)), 4, 8, fileNames((pDs==1) & (abs(TI)>0.1)));
    xlim([0,1]);
    refline(1);
    refline(0,0.5);
    reflinexy(0.5,1);
    xlabel('Grand CP');
    ylabel('Idisp Effect');
    
    figure(9669), clf, hold on,
    scatterhist(GrandCP(GrandCP>0), IdBiasROC1(GrandCP>0), ...
        'NBins' , [20, 20], ...
        'Location', 'SouthEast', ...
        'Direction' , 'out');
    xlim([0,1]);
    ylim([0,1]);
    refline(1);
    refline(0,0.5);
    reflinexy(0.5,1);
    xlabel('Grand CP');
    ylabel('Idisp Effect');
    [r, p] = corr(GrandCP(GrandCP>0)', IdBiasROC1(GrandCP>0)');
    text(0.05,0.9, ['Correlation  r:  ' num2str(r) '   p: ' num2str(p)]);
    
end

%% DPI

if (strcmpi(FileType, 'DPI'))
    figure(6785), clf, clickscatter(TIcyldx, PI, 6, 7, fileNames);
    xlabel('Tuning Index');
    ylabel('Pursuit Index : (Fr(pref) - Fr(null) / Fr(pref) + Fr(null))');
    
    figure(6755), clf, hold on
    clickscatter(TI, PIPrefdx, 7, 8, fileNames);
    clickscatter(TI, PINulldx, 6, 8, fileNames);
    xlabel('Tuning Index');
    ylabel('Pursuit Index : (Fr(pref) - Fr(null) / Fr(pref) + Fr(null))');
    
    figure(6665), clf, hold on
    clickscatter(TI, PINew, 8, 8, fileNames);
    clickscatter(TI, PIPrefdx, 7, 8, fileNames);
    clickscatter(TI, PINulldx, 6, 8, fileNames);
    xlabel('Tuning Index');
    ylabel('Pursuit Index : (Fr(pref) - Fr(null) / Fr(pref) + Fr(null))');
    
    figure(1982), clf, hold on,
    clickscatter(abs(TI), abs(PINew), 8, 8, fileNames);
    clickscatter(abs(TI(abs(TI)>0.1)), abs(PINew(abs(TI)>0.1)), 6, 8, fileNames(abs(TI)>0.1));
    xlabel('ABS( Tuning Index )');
    ylabel('ABS( Pursuit Index) : ABS((Fr(pref) - Fr(null) / Fr(pref) + Fr(null)) )');
    
    
    %%
    figure(1357), clf, hold on,
    %     bar([mean(abs(PIPrefdx(~isnan(PIPrefdx)))), mean(abs(PIZerodx(~isnan(PIZerodx)))), mean(abs(PINulldx(~isnan(PINulldx)))), mean(abs(PI(~isnan(PI))))]);
    %     errorbar([mean(abs(PIPrefdx(~isnan(PIPrefdx)))), mean(abs(PIZerodx(~isnan(PIZerodx)))), mean(abs(PINulldx(~isnan(PINulldx)))), mean(abs(PI(~isnan(PI))))], ...
    %              [std(abs(PIPrefdx(~isnan(PIPrefdx))))./sqrt(sum(~isnan(PIPrefdx))), std(abs(PIZerodx(~isnan(PIZerodx))))./sqrt(sum(~isnan(PIZerodx))), std(abs(PINulldx(~isnan(PINulldx))))./sqrt(sum(~isnan(PINulldx))), std(abs(PI(~isnan(PI))))./sqrt(sum(~isnan(PI)))]);
    %
    %     [a, b, c] = anova1([abs(PI)', abs(PIPrefdx)', abs(PINulldx)', abs(PIZerodx)']);
    PiP = PIPrefdx .* (TI > 0) - PIPrefdx .* (TI < 0);
    PiN = PINulldx .* (TI > 0) - PINulldx .* (TI < 0);
    PiZ = PIZerodx .* (TI > 0) - PIZerodx .* (TI < 0);
    PiA = PINew    .* (TI > 0) - PINew    .* (TI < 0);
    
    bar([mean(abs(PiP(~isnan(PiP)))), ...
        mean(abs(PiZ(~isnan(PiZ)))), ...
        mean(abs(PiN(~isnan(PiN)))), ...
        mean(abs(PiA(~isnan(PiA)))), ...
        mean(abs(PiA(~isnan(PiA) & abs(TI)<0.1)))]);
    errorbar([mean(abs(PiP(~isnan(PiP)))), ...
        mean(abs(PiZ(~isnan(PiZ)))), ...
        mean(abs(PiN(~isnan(PiN)))), ...
        mean(abs(PiA(~isnan(PiA)))), ...
        mean(abs(PiA(~isnan(PiA) & abs(TI)<0.1)))], ...
        [std(abs(PiP(~isnan(PiP))))./sqrt(sum(~isnan(PiP))), ...
        std(abs(PiZ(~isnan(PiZ))))./sqrt(sum(~isnan(PiZ))), ...
        std(abs(PiN(~isnan(PiN))))./sqrt(sum(~isnan(PiN))), ...
        std(abs(PiA(~isnan(PiA))))./sqrt(sum(~isnan(PiA))), ...
        std(abs(PiA(~isnan(PiA) & abs(TI)<0.1)))./sqrt(sum(~isnan(PiA) & abs(TI)<0.1))]);
    
    [a, b, c] = anova1([abs(PiP)', abs(PiZ)', abs(PiN)', abs(PiA)']);
    multcompare(c);
    
    %%
    PISliding = [];
    TIs = sort(TI);
    for slide=1:length(TI)-10
        PISliding(slide) = mean(PiA((TI>=TIs(slide)) & (TI<TIs(slide+10))));
    end
    
    figure(1999), clf, hold on
    plot(TIs(1:end-10),PISliding);
    refline(0,0);
    reflinexy(0,1);
    
end

%% crosstalk
if (exist('bdCrossTalk5')==1)
    if ~strcmpi(FileType, 'DPI') && ~strcmpi(FileType, 'JPI')
        figure(1819), hist(bdCrossTalk5(abs(orS)==90))
        hold on , hist(bdCrossTalk5((orS==0) | (abs(orS)==180)))
    end
end

%% or
if ~strcmpi(FileType, 'DPI') && ~strcmpi(FileType, 'JPI')
    bb = [0 45 90 135 ];
    for i = 1: length(orS),
        switch round(orS(i))
            case {-90, 80, 110, 115}, orSm(i) = 90;
            case 180, orSm(i) = 0;
            case {-135, -160, 30}, orSm(i) = 45;
            case {-30, -45, 120}, orSm(i) = 135;
            otherwise orSm(i) = orS(i);
        end,
    end
    g = zeros(length(orSm),1);
    for i = 1: length(bb), g(orSm==bb(i)) = i; end
    
    m = []; e= [];
    for i = 1:length(bb), m(i) = mean(IdBiasROC1(round(orSm) == round(bb(i)))); end
    for i = 1:length(bb), e(i) = std(IdBiasROC1(round(orSm) == round(bb(i)))); end
    
    figure(856), clf
    errorbar(bb, m, e)
    figure(678), scatter(orSm, IdBiasROC1, 'filled')
    figure(698), scatter(orSm, BiaEff, 'filled')
    
    [p,table,stats] = anova1(IdBiasROC1(abs(TI)>0.1), g(abs(TI)>0.1))
    [p,table,stats] = anova1(BiaEff(abs(TI)>0.1), g(abs(TI)>0.1))
    
end


%% Graphics

switch(upper(FileType))
    case 'JPI'
        figure(816), clf, clickscatter(TI, dpJPIRDS1, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
        figure(826), clf, clickscatter(TI, dpJPIRDS2, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
        figure(836), clf, clickscatter(TI, dpJPIRDS3, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
        figure(846), clf, clickscatter(TI, dpJPIRDS4, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
        
    case 'DPI'
        if strcmpi(StimulusType, 'cylinder')
            
            for i = 1: length(dprimes),
                if(~isempty(dprimes{i}))
                    if(~isnan(dprimes{i}(1)))
                        for j = 1: length(dprimes{i})
                            dps(i,j) = dprimes{i}(j);
                        end
                    end
                end
            end
            
            figure(326), clf,
            clickscatter(TI, dps(:,4)', 6 + pDs, 8, fileNames);
            refline(0, 0); reflinexy(0,100);
            
            % d prime shift
            figure(336), clf,
            clickscatter(TI, (dps(:,4) - ((dps(:,8) + dps(:,12))./2))', 8, 7, fileNames);
            refline(0, 0); reflinexy(0,100);
            %  normalized d prime shift
            figure(346), clf,
            scatter(TI, ((dps(:,8) - dps(:,12))./ (dps(:,8) + dps(:,12)))', DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
            refline(0, 0); reflinexy(0,100);
            
            figure(356), clf,
            clickscatter(TI, dpdxless, 6 + pDs, 8, fileNames);
            refline(0, 0); reflinexy(0,100);
            
            figure(366), clf,
            clickscatter(TI, dpdx, 6 + pDs, 8, fileNames);
            refline(0, 0); reflinexy(0,100);
            
            % average d prime
            figure(386), clf,
            clickscatter(TI, (dps(:,4) + dps(:,8) + dps(:,12))'./3, 8, 7, fileNames);
            refline(0, 0); reflinexy(0,100);
            
        else % rds
            figure(616), clf, clickscatter(TI, dpRDS1, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
            figure(626), clf, clickscatter(TI, dpRDS2, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
            figure(636), clf, clickscatter(TI, dpRDS3, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
            figure(646), clf, clickscatter(TI, dpRDS4, 6 + (TI>0), 8, fileNames); refline(0, 0); reflinexy(0,100);
        end
    case 'BDID'
        figure(1123), clf, hist(IdBiasROC2);
        figure(1145), clf, clickscatter(TI, IdBiasROC2, 1, 7, fileNames); refline(0, 0.5);
    case 'TWO'
        figure(190), clf, clickscatter(IdBiasROC1, FlipROC, 1+(BiaEff>BiaEffFlip), 7, fileNames); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        ylabel('Flip ROC');
        xlabel('Main ROC');
        ylim([0.2 0.8]);
        xlim([0.2 0.8]);
        refline(0, 0.5);
        reflinexy(0.5, 1);
        
        
        figure(291), clf, clickscatter(IdBiasROC1(GrandCPFlip>0), GrandCPFlip(GrandCPFlip>0), 1+(BiaEff>BiaEffFlip), 7, fileNames); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        ylabel('GrandCP for flip');
        xlabel('Main ROC');
        ylim([0.1 1.]);
        xlim([0.1 1.]);
        refline(0, 0.5);
        reflinexy(0.5, 1);
        
        figure(298), clf, clickscatter(FlipROC(GrandCPFlip>0), GrandCPFlip(GrandCPFlip>0), 1+(BiaEff>BiaEffFlip), 7, fileNames); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        ylabel('GrandCP for flip');
        xlabel('Flip ROC');
        ylim([0.1 1.]);
        xlim([0.1 1.]);
        refline(0, 0.5);
        reflinexy(0.5, 1);
        
        
        
        figure(842), clf,
        clickscatter(RFproximity, FlipROC, 1+(BiaEff>BiaEffFlip), 7, fileNames); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        refline(0, 0.5);
        reflinexy(0.5, 1);
        ylabel('Flip ROC');
        xlabel('RF proximity');
        ylim([0.3 0.8]);
        %xlim([0.2 0.9])
        
    otherwise
        figure(17), clf, clickscatter(TI, IdBiasROC1, pDs, 7, fileNames) %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        refline(0, 0.5);
        %figure, scatter(abs(TI), IdBiasROC1, [], reshape(([IdColor{:}]), 3,length(IdBiasROC1))', 'filled');
        
        sum(BiaEff>0) / length(BiaEff) % percentage inclusion for Correctly biased thing
        figure(23), clf, clickscatter(TI(BiaEff>0), IdBiasROC1(BiaEff>0),1, 7, fileNames); %, DotSizes(BiaEff>0), reshape(([IdColor{BiaEff>0}]), 3, sum(BiaEff>0))', 'filled');
        refline(0, 0.5);
        reflinexy(0,1);
        
end



%%
figure(22222), hist(GrandCP(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0), [0.05:0.1:0.95])
figure(22222), hold on, hist(GrandCP(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0 & GrandCPSig<5), [0.05:0.1:0.95])
figure(11111), hist(IdBiasROC1(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0), [0.05:0.1:0.95])
figure(11111), hold on, hist(IdBiasROC1(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0 & IdBiasROCSig<5), [0.05:0.1:0.95])
%%
ax = [0.0:0.05:1];
awhite = hist(GrandCP(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0), ax);
ablack = hist(GrandCP(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0 & GrandCPSig<5), ax);
figure(220022), clf, hold on,
bar(ax, awhite+ ablack, 'r');
bar(ax, ablack, 'k');
xlim([0 1]);
%%
bwhite = hist(IdBiasROC1(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0), ax);
bblack = hist(IdBiasROC1(abs(TI)>0.05 & IdBiasROC1>0 & GrandCP>0 & IdBiasROCSig<5), ax);
figure(110011), clf, hold on,
bar(ax, bwhite + bblack, 'r');
bar(ax, bblack, 'k');
xlim([0 1]);


%%
trialcountstable = zeros(2,2);
zscoredspikecountstable = zeros(2,2);
cellcounter = 0;
for i = 1:length(ContignGTables)
    if shart(i)
         if (sum(isnan(ContignGTables{i}(:)))==0)
            trialcountstable = trialcountstable + ContignGTables{i}(:,:,2);
            zscoredspikecountstable = zscoredspikecountstable + ContignGTables{i}(:,:,1);
            cellcounter = cellcounter + 1;
         end
    end
end

disp(['cells in this table: ', num2str(cellcounter)]);

%model
%
% ----------------------------------------------------------------------
%       |  p choice                     | null choice                  |
% ------|-------------------------------|------------------------------|
% p bd  | ((b+t)(n-m1)+2tm1) / (n+m1)   |  -b-t                        |
% ------|-------------------------------|------------------------------|
% n bd  | b + t                         | ((-b-t)(n-m2)-2tm2) / (n+m2) |
% ----------------------------------------------------------------------

% b+t -> PN
PN = somesortofnormalizedorweightedzscoredspikecountstable(2,1);
% -b-t -> NP
NP = somesortofnormalizedorweightedzscoredspikecountstable(1,2);
%  ((b+t)(n-m1)+2tm1) / (n+m1)



%%
load ~/Downloads/twoall.mat
figure(44444), hist(DATA.Data.ally(DATA.Data.selected==1), [0.25:0.1:0.75]); xlim([0.2 0.8])
figure(44444), hold on, hist(DATA.Data.ally(DATA.Data.selected==1 & (DATA.Data.rocpval(:,2)<0.05)'), [0.25:0.1:0.75]); xlim([0.2 0.8])
figure(33333), hist(DATA.Data.allx(DATA.Data.selected==1), [0.25:0.1:0.75]); xlim([0.2 0.8])
figure(33333), hold on, hist(DATA.Data.allx(DATA.Data.selected==1 & (DATA.Data.rocpval(:,1)<0.05)'), [0.25:0.1:0.75]); xlim([0.2 0.8])



%% ROC distribution analysis
% is our ROCs come from a single distribution or two or more groups of neurons

ROCvarTest(ROCpairs1);
if strcmp(FileType, 'TWO')
    ROCvarTest(ROCpairsFlip);
    % else
    %     if strcmp(FileType, 'BDID')
    %         ROCvarTest(ROCpairsFlip);
    %     end
end


%% Compaare with Bruce's
% figure, scatter(IdBiasROC1, flipud(data(:,2)), DotSizes, reshape(([IdColor{:}]), 3,63)', 'filled');
%
% data = flipud(data);
% figure, scatter(TI, data(:,1), DotSizes, reshape(([IdColor{:}]), 3, 63)', 'filled')

