% This file has been cross check with B. and is absolutely bug free as
% of 10/16/2009 for the AllABDNeurons.txt for the ABD and cylinder
% obviously.

clear;
clc

DataPath = GetDataPath();

FinishTime = 0;

% ABD
% load /Users/ali/DropBox/Projects/BCode/AllABDNeurons.mat
% AllNeurons = AllABDNeurons;
%  % AllNeurons = AllABDNeurons(57:end); disp ( ' O N L Y   I C A R U S ! ! ! !   B E W A R E ! ! !');
% clear AllABDNeurons;
% FileType = 'ABD';
% StimulusType = 'cylinder';
% StartTime  = 10000;  %10000; % 6500;
% FinishTime = 20000;

% % DID
load /Users/ali/DropBox/Projects/BCode/AllDIDNeurons.mat
AllNeurons = AllDIDNeurons;
clear AllDIDNeurons
FileType = 'DID';
StimulusType = 'cylinder';
StartTime  = 10000; %10000; % 6500; 
FinishTime = 20000;

% % TWO
% load('/Users/ali/DropBox/Projects/BCode/AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWONeurons;
% FileType = 'TWO';
% StimulusType = 'cylinder';
% StartTime  = 500; %10000; %5500; 
% FinishTime = 20000;

% % DPI
% load('/Users/ali/DropBox/Projects/BCode/AllPursuitNeurons.mat');
% AllNeurons = AllPursuitNeurons;
% clear AllPursuitNeurons;
% FileType = 'DPI';
% StimulusType = 'cylinder';
% StartTime  = 500;% 10000; %500; %10000; 
% FinishTime = 20000;


% % % BDID
% load('/Users/ali/DropBox/Projects/BCode/AllBDIDNeuronsALL.mat');
% AllNeurons = AllBDIDNeuronsALL;
% clear AllBDIDNeurons;
% FileType = 'BDID';
% StimulusType = 'cylinder';
% StartTime  = 10000;
% FinishTime = 20000;
% 


%par
for iN= [1 :length(AllNeurons)], 
    if iN == 38
        debug = 1;
    end
    IdColor{iN} = [1 0.5 0.1];
    DotSizes(iN) = 100;
    toolowFR = 0;
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
    
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    fileNames{iN} = filename; 
    
    pD      = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    %if TI(iN) > 0, pD = 1; else pD = 2; end
    pDs(iN) = pD;

    if (pD == -1) 
        disp('Ohoooy! pD? what are you D O I N G ! ? ');
    end
    
    %conditions = GetConditions(Expt, FileType, pD);

    conditions = logical([]);
    
    if strcmpi(FileType, 'DPI')
         if isfield(Expt.Trials,'dfx')
            deltafxy = [Expt.Trials(:).dfx] - [Expt.Trials(:).fx];
        else
            deltafxy = [Expt.Trials(:).dfy] - [Expt.Trials(:).fy];
        end
        tempdelta = fix(deltafxy);
        for td = 1: length(tempdelta)
            if tempdelta(td) == 0
                tempdelta(td) = 0.5 * sign(deltafxy(td));
            end
        end
        deltafxy = tempdelta; % fix(deltafxy);  %  0.2 * fix(5*deltafxy); %0.1 * round(10*deltafxy);


        Speeds = unique(abs(deltafxy));
        %disp(length(Speeds));
        %disp(Speeds);
       if(Expt.Stimvals.or>0)
            if(pD==2)
                conditions(1,:) = (deltafxy>0 & [Expt.Trials(:).dx]==0);
                conditions(2,:) = (deltafxy<0 & [Expt.Trials(:).dx]==0);
                conditions(3,:) = (deltafxy>0 & [Expt.Trials(:).dx]>0);
                conditions(4,:) = (deltafxy<0 & [Expt.Trials(:).dx]>0);
                conditions(5,:) = (deltafxy<0 & [Expt.Trials(:).dx]<0);
                conditions(6,:) = (deltafxy>0 & [Expt.Trials(:).dx]<0);
            else
                conditions(1,:) = (deltafxy<0 & [Expt.Trials(:).dx]==0);
                conditions(2,:) = (deltafxy>0 & [Expt.Trials(:).dx]==0);
                conditions(3,:) = (deltafxy<0 & [Expt.Trials(:).dx]>0);
                conditions(4,:) = (deltafxy>0 & [Expt.Trials(:).dx]>0);
                conditions(5,:) = (deltafxy>0 & [Expt.Trials(:).dx]<0);
                conditions(6,:) = (deltafxy<0 & [Expt.Trials(:).dx]<0);
            end
       else
            if(pD==1)
                conditions(1,:) = (deltafxy>0 & [Expt.Trials(:).dx]==0);
                conditions(2,:) = (deltafxy<0 & [Expt.Trials(:).dx]==0);
                conditions(3,:) = (deltafxy>0 & [Expt.Trials(:).dx]>0);
                conditions(4,:) = (deltafxy<0 & [Expt.Trials(:).dx]>0);
                conditions(5,:) = (deltafxy<0 & [Expt.Trials(:).dx]<0);
                conditions(6,:) = (deltafxy>0 & [Expt.Trials(:).dx]<0);
            else
                conditions(1,:) = (deltafxy<0 & [Expt.Trials(:).dx]==0);
                conditions(2,:) = (deltafxy>0 & [Expt.Trials(:).dx]==0);
                conditions(3,:) = (deltafxy<0 & [Expt.Trials(:).dx]>0);
                conditions(4,:) = (deltafxy>0 & [Expt.Trials(:).dx]>0);
                conditions(5,:) = (deltafxy>0 & [Expt.Trials(:).dx]<0);
                conditions(6,:) = (deltafxy<0 & [Expt.Trials(:).dx]<0);
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
    

        % 10/14/09 we added the [Expt.Trials(:).RespDir]~=0 to all conditions
        % (firstly to match is with Bruces secondly to exclude those in which
        % monkey has not taken a side. 
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
            end
     else if strcmp(FileType, 'BDID')
             conditions = GetConditions(Expt, FileType, pD);
%              if(pD == 2)
%                     conditions(1,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(2,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(3,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                     conditions(4,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                     conditions(5,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                     conditions(6,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                     conditions(7,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0; 
%                     conditions(8,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0; 
%                     conditions(9,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]~=0; 
%                     conditions(10,:)= [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]~=0; 
%                 else
%                     conditions(1,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(2,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(3,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                     conditions(4,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                     conditions(5,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                     conditions(6,:) = [Expt.Trials(:).bd]<0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                     conditions(7,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(8,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0;
%                     conditions(9,:) = [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]~=0; 
%                     conditions(10,:)= [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]~=0; 
%              end
         else
             conditions = GetConditions(Expt, FileType, pD);
%             if(pD == 2) % DID , ...
%                 conditions(1,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
%                 conditions(2,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
%                 conditions(7,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                 conditions(8,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                 conditions(9,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                 conditions(10,:)= [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                 conditions(11,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                 conditions(12,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%             else 
%                 conditions(1,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
%                 conditions(2,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]~=0;
%                 conditions(7,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                 conditions(8,:) = [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                 conditions(9,:) = [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%                 conditions(10,:)= [Expt.Trials(:).Id]>0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                 conditions(11,:)= [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative;
%                 conditions(12,:)= [Expt.Trials(:).Id]<0 & [Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive;
%             end
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
        SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
    end
    
    if (strcmpi(FileType, 'DPI'))
        dprs = []; dprsdxp=[]; dprsdxn=[];
        for ss = 1: length(Speeds),
            cs = (abs(deltafxy) == Speeds(ss));
            dprs(ss) = dPrime(SpikeCounts(conditions(1,:) & cs), SpikeCounts(conditions(2,:) & cs));
            dprsdxp(ss) = dPrime(SpikeCounts(conditions(3,:) & cs), SpikeCounts(conditions(4,:) & cs));
            dprsdxn(ss) = dPrime(SpikeCounts(conditions(6,:) & cs), SpikeCounts(conditions(5,:) & cs));
        end
        dprs(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
        dprsdxp(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(4,:)));
        dprsdxn(1+ length(Speeds)) = dPrime(SpikeCounts(conditions(6,:)), SpikeCounts(conditions(5,:)));

        dprimes{iN} = [dprs dprsdxp dprsdxn];

        NS{iN} = {SpikeCounts, conditions, deltafxy};
 
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
    else % DID, ABD, TWO ...
        IdBiasROC1(iN) = ROCAUC(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:)));
        [rr, pp] = ROCAUCSignificance(SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:))); 
        IdBiasROCSig(iN) = pp;
        ROCpairs1{iN} = {SpikeCounts(conditions(1,:)), SpikeCounts(conditions(2,:))};
    end
    switch FileType
        case {'TWO', 'BDID'}
            if Expt.Stimvals.bo == Expt.Stimvals.or
                Manip(iN) = 1;
            else
                Manip(iN) = 2;
            end
        otherwise
            CP(iN)      = ROCAUC(SpikeCounts(conditions(7,:)), SpikeCounts(conditions(8,:)));
            CPPref(iN)  = ROCAUC(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(5,:))); %CPPref(iN)  = ROCAUC(SpikeCounts(conditions(9,:)), SpikeCounts(conditions(10,:)));
            CPNull(iN)  = ROCAUC(SpikeCounts(conditions(4,:)), SpikeCounts(conditions(6,:))); %CPNull(iN)  = ROCAUC(SpikeCounts(conditions(11,:)), SpikeCounts(conditions(12,:)));
            %CP1(iN)     = ROCAUC(SpikeCounts(conditions(3,:)), SpikeCounts(conditions(5,:)));
            %CP2(iN)     = ROCAUC(SpikeCounts(conditions(6,:)), SpikeCounts(conditions(4,:)));
            %CPm(iN)     = mean([CP1(iN) CP2(iN)]);
            %CPeffect(iN)= CP1(iN) - CP2(iN);
    end
    if strcmp(FileType,'TWO')
        fa = ROCAUC(SpikeCounts(conditions(8,:)), SpikeCounts(conditions(7,:)));
        fb = ROCAUC(SpikeCounts(conditions(10,:)), SpikeCounts(conditions(11,:)));
        if abs(fa * fb ) ~= (fa * fb)
            disp('Fishy - - - - - - - - - - - 9 9 9 9 9 9 8  8 8 8( ( ((  ** * * * ');
            fa = 0; fb = 0;
        end
        FlipROC(iN) = 0.5 * ( fa + fb );
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
        else
            BiaEff(iN) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0) ...
               - 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
        end
        end
    end    

    if IdBiasROC1(iN) > 0.7
        disp( IdBiasROC1(iN));
        debug = 1;
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
end


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
%% Graphics

if(strcmpi(FileType, 'DPi'))
    for i = 1: length(dprimes),
        if(~isempty(dprimes{i}))
            if(~isnan(dprimes{i}(1)))
                for j = 1: length(dprimes{i})
                    dps(i,j) = dprimes{i}(j);
                end
            end
        end
    end
    
    figure, clf,
    brucescatter(TI, dps(:,4)') %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
    refline(0, 0); reflinexy(0,100);

    % d prime shift
    figure, clf,
    brucescatter(TI, (dps(:,4) - ((dps(:,8) + dps(:,12))./2))'); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
    refline(0, 0); reflinexy(0,100);
    %  normalized d prime shift
    figure, clf,
    scatter(TI, ((dps(:,8) - dps(:,12))./ (dps(:,8) + dps(:,12)))', DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
    refline(0, 0); reflinexy(0,100);
   
else
    if(strcmpi(FileType, 'BDID'))
        figure(1123), clf, hist(IdBiasROC2);
        figure(1145), clf, brucescatter(TI, IdBiasROC2, 1, 7, fileNames); refline(0, 0.5);
    else
    figure(17), clf, brucescatter(TI, IdBiasROC1, pDs, 7, fileNames) %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
    refline(0, 0.5);
    %figure, scatter(abs(TI), IdBiasROC1, [], reshape(([IdColor{:}]), 3,length(IdBiasROC1))', 'filled');
    if (strcmp(FileType, 'TWO'))
        figure(19), clf, brucescatter(FlipROC, IdBiasROC1, 1+(BiaEff>BiaEffFlip), 7, fileNames); %, DotSizes, reshape(([IdColor{:}]), 3,length(IdColor))', 'filled');
        refline(0, 0.5);
    end

    sum(BiaEff>0) / length(BiaEff) % percentage inclusion for Correctly biased thing
    figure(23), clf, brucescatter(TI(BiaEff>0), IdBiasROC1(BiaEff>0),1, 7, fileNames); %, DotSizes(BiaEff>0), reshape(([IdColor{BiaEff>0}]), 3, sum(BiaEff>0))', 'filled');
    refline(0, 0.5);
    reflinexy(0,1);
    end
end
%% Compaare with Bruce's
% figure, scatter(IdBiasROC1, flipud(data(:,2)), DotSizes, reshape(([IdColor{:}]), 3,63)', 'filled');
% 
% data = flipud(data);
% figure, scatter(TI, data(:,1), DotSizes, reshape(([IdColor{:}]), 3, 63)', 'filled')

