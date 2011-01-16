clear;
clc;

ShowSingleCellSDFs = 0; % 0 or 1

% % ABD
% load('/Users/ali/DropBox/Projects/BCode/AllABDNeurons.mat');
% AllNeurons = AllABDNeurons;
% clear AllABDNeurons;
% FileType = 'ABD';
% StimulusType = 'cylinder';

% % % DID
% load('/Users/ali/DropBox/Projects/BCode/AllDIDNeurons.mat');
% AllNeurons = AllDIDNeurons;
% clear AllDIDNeurons;
% FileType = 'DID';
% StimulusType = 'cylinder';


% % % TWO
% load('/Users/ali/DropBox/Projects/BCode/AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWONeurons;
% FileType = 'TWO';
% StimulusType = 'cylinder';

% % BDID
load('/Users/ali/DropBox/Projects/BCode/AllBDIDNeuronsALL.mat');
AllNeurons = AllBDIDNeuronsALL;
clear AllBDIDNeuronsALL;
FileType = 'BDID';
StimulusType = 'cylinder';
% AllNeurons = AllNeurons(1:37); 
% %disp('= = = = = =  JUST LOOKING AT Adrian s  Cells = = = = = =');
%  AllNeurons = AllNeurons(38:51); 
%  disp('= = = = = =  Ignoring  Adrian s  Cells = = = = = =');

% % % DIDB
% load('/Users/ali/DropBox/Projects/BCode/AllDIDBNeurons.mat');
% AllNeurons = AllDIDBNeurons;
% clear AllDIDBNeurons;
% FileType = 'DIDB';
% StimulusType = 'cylinder';

% % DRID
% load('/Users/ali/DropBox/Projects/BCode/AllDRIDNeurons.mat');
% FileType = 'DRID';
% AllNeurons = AllDRIDNeurons;
% % SRID
% load('/Users/ali/DropBox/Projects/BCode/AllSRIDNeurons.mat');
% FileType = 'SRID';
% AllNeurons = AllSRIDNeurons;
% clear AllDRIDNeurons AllSRIDNeurons;
% StimulusType = 'rds';

%Prep
DataPath = '/bgc/data/';
BinSize = 50;%50;
SmoothingBinSize = 1;%50
SmthKernel = gausswin(SmoothingBinSize);

filenamesforbruce = {};
TI=[];

% THIS IS THE SELETION OF DRID dx tunied experiemtns
%AllNeurons = AllNeurons([3,4,6,7, 11, 13, 15, 20, 21, 22,23]);
%par
parfor iN= 1:length(AllNeurons) %[length(AllNeurons):-1:1], 1:length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));

    p = []; c = []; eb = []; jnk1 = []; jnk2 = []; xC =[]; c1 =[]; ebX =[];
    
    switch FileType 
        case {'ABD', 'DID', 'DRID', 'SRID', 'DIDB'} 
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
            %if (abs(TI(iN))>0.1)
                %[xC ,c1, ebX] = PlotSTAutoCorr(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, ShowSingleCellSDFs);
                [p, c, eb, jnk1, jnk2] = PlotPSTHIdisp(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
                %[c, eb] = PlotISISdf(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, 1);
                slp{iN} = PsychSlop(AllNeurons(iN), StimulusType, FileType, 'dx'); 
                %disp(['Ti; ', num2str(TI(iN)), ' - Slop: ' , num2str(slp(iN).fit(2))]);
            %end
        case 'BDID'
            [p, c, eb] = PlotPSTHIdisp(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
            slp{iN} = PsychSlop(AllNeurons(iN), StimulusType, FileType, 'bd'); 
            %disp(['Ti: ', num2str(TI(iN)), ' - Slop: ' , num2str(slp(iN).fit(2))]);
        case 'TWO'
            [p, c, eb] = PlotPSTHIdisp(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
            %[p, c, eb] = PlotPSTHTWO(NeuronNumber, FileType, StimulusType, BinSize, 0);
            %TI(iN) = TuningIndex(NeuronNumber, ClusterName, StimulusType, FileType, 'dx');
          
    end
    
    %eb = convn(eb, SmthKernel')./ sum(SmthKernel);
    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    if (strcmp(FileType, 'DRID') || strcmp(FileType, 'SRID'))
        cellValid  = 1; vScore = 1;
    else
        [cellValid , vScore] = CellValidity([], FileType, pD, MonkeyName, NeuronNumber, ClusterName, FileType);
    end
    PSTHs{iN} = {p , c, eb, cellValid, vScore, pD, jnk1, jnk2, xC, ebX};
    filenamesforbruce{iN} = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
end

%%
validCells = []; for i = 1: length(PSTHs), if~isempty(PSTHs{i}), validCells(i) = [PSTHs{i}{4}]; end; end
AllConditions = [];

for i = 1 : length(AllNeurons),
    if (~isempty(PSTHs{i}) && ~isempty(PSTHs{i}{1}))
        if ~isempty(PSTHs{i}{7}), rocs(i,:) = PSTHs{i}{7}; end
        if ~isempty(PSTHs{i}{8}), DandC(i,:,1:2100) = PSTHs{i}{8}(:,1:2100); end
        %xCs(i,:,:) = PSTHs{i}{9};
        ebX = (PSTHs{i}{10});
        p = ([PSTHs{i}{1}]);
        if (isempty(AllConditions)) AllConditions = zeros(length(AllNeurons), length(sum([PSTHs{i}{2}],2))); end
        AllConditions(i,:) = sum([PSTHs{i}{2}],2);
        eb = ([PSTHs{i}{3}]);
        pD(i) = ([PSTHs{i}{6}]);
        
        cntr = 0;
        for ebi = 1:size(eb,1)
            cntr = cntr + 1;
            if(sum(eb(ebi,:))>0)
                tPSTH(i, cntr,1:size(eb(cntr,:),2)) = eb(ebi,:);
                if(ebi<=size(ebX,1))
                    txCor(i, cntr,1:size(ebX(cntr,:),2))= ebX(ebi,:);
                end
            end
        end
    end
end

%%
switch FileType
    case {'DID', 'ABD', 'DIDB'}
        if ~isempty(slp)
           % vs = reshape([slp.fit], 2, length(AllNeurons));
           for i = 1:length(slp), vs(i) = slp{i}.fit(2); end
        end
        criteria = ((abs(TI)>0.1)); 
    case 'BDID'
        if ~isempty(slp)
            for i = 1:length(slp), vs(i) = slp{i}.fit(2); end
            %vs = reshape([slp.fit], 2, length(AllNeurons));
        end
        %criteria = ((abs(TI)>0.1) & (vs(2,:) >0) & (vs(2,:)<0.01));
        criteria = ((abs(TI)>0.1) & (vs >0) & (vs<0.01));
    case 'TWO'
        criteria = (TI>0.1 | TI<-0.1);% & validCells;
end

%% Graphics 

%PopPSTH = mean(nPSTH);
%PopPSTH = mean(tPSTH(~isnan(mean(mean(tPSTH,2),3))' ,:,:));
if strcmp(FileType,'TWO')
    %PopPSTH = mean(tPSTH(~isnan(mean(mean(tPSTH,2),3))' & (TI>0.1 | TI<-0.1) & validCells,:,:));
    %PopPSTH = mean(tPSTH(~isnan(mean(mean(tPSTH,2),3))' & (TI>0.1 | TI<-0.1),:,:));
    PopPSTH = mean(tPSTH(criteria,:,:));
else
%    PopPSTH = mean(tPSTH(~isnan(mean(mean(tPSTH,2),3))' & (TI>0.1 | TI<-0.1),:,:));
    %PopPSTH = mean(tPSTH(~isnan(mean(mean(tPSTH,2),3))' & (TI>0.05 | TI<-0.05),:,:));
    PopPSTH = squeeze(mean(tPSTH(criteria,:,:)));
    PopPSTHse = std(tPSTH(criteria,:,:))/sqrt(size(tPSTH,1));
    if exist('txCor')
        PopxCor = squeeze(mean(txCor(criteria,:,:)));
    end
end



%%
w_mode = 'trial_count';%'pref_trial'; %'pref'; %'trial_count';
switch w_mode
    case 'trial_count'
        for wi =1: size(tPSTH,1)
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = AllConditions(wi,wc) * squeeze(tPSTH(wi,wc,:));
            end
        end
        WeightedPopPSTH = squeeze(sum(WtPSTH(criteria,:,:)));
        for i = 1: size(WeightedPopPSTH,2)
            WeightedPopPSTH(:,i) = WeightedPopPSTH(:,i) ./ sum(AllConditions(criteria,:))';
        end
    case 'pref'
        for wi =1: size(tPSTH,1)
            prefMean = mean(squeeze(tPSTH(wi, 1, :)));
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = squeeze(tPSTH(wi,wc,:))/prefMean;
            end
        end
        WeightedPopPSTH = squeeze(mean(WtPSTH(criteria,:,:)));
    case 'pref_trial'
        for wi =1: size(tPSTH,1)
            prefMean = mean(squeeze(tPSTH(wi, 1, :)));
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = squeeze(tPSTH(wi,wc,:))/prefMean;
            end
        end
        for wi =1: size(tPSTH,1)
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = AllConditions(wi,wc) * squeeze(WtPSTH(wi,wc,:));
            end
        end
        WeightedPopPSTH = squeeze(sum(WtPSTH(criteria,:,:)));
        for i = 1: size(WeightedPopPSTH,2)
            WeightedPopPSTH(:,i) = WeightedPopPSTH(:,i) ./ sum(AllConditions(criteria,:))';
        end

end



%%
switch FileType
    case 'SRID'
        figure(142);
    case 'DRID'
        figure(143);
    otherwise
        figure(125);
end
clf, hold on,  
h = plot(squeeze(PopPSTH)');
%h = plot(squeeze(WeightedPopPSTH)');
if strfind(FileType, 'RID')
    plot(PopPSTH(1,:) - PopPSTH(2,:), 'm', 'LineWidth', 2);
    plot(PopPSTH(3,:) - PopPSTH(4,:), 'k', 'LineWidth', 2);
    refline(0,0);
end

set(h, 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([100 2200]);
xtl = [0, 50, 500, 1000, 1500, 2000];
set(gca, 'XTick', xtl+200-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);

%errorbar(squeeze(PopPSTH)', squeeze(PopPSTHse)');

sum(criteria)

% a = squeeze(PopPSTH)';
% %hold on, plot(CumulativeDifference(a(:,1), a(:,2))/100)
% hold on, 
% 
% if strcmp(FileType, 'TWO')
%     plot(a(:,3) - a(:,4), 'r');
% else
%     plot(a(:,1) - a(:,2), 'r:');
%     plot(a(:,3) - a(:,4), 'r');
% end

%% Psych Effect vs Neuronal Id Effect

switch FileType
    case {'DIDB', 'DID'}
        for i = 1: length(AllNeurons), 
            if (pD(i) == 2), 
                PsychEff(i) = ((AllConditions(i,3) / (AllConditions(i,3) + AllConditions(i,5))) - (AllConditions(i,6) / (AllConditions(i,4) + AllConditions(i,6)))); 
                NeuIdEff(i) = ((mean(tPSTH(i,3,1200:2200),3) / (mean(tPSTH(i,3,1200:2200),3) + mean(tPSTH(i,5,1200:2200),3))) - (mean(tPSTH(i,6,1200:2200),3) / (mean(tPSTH(i,4,1200:2200),3) + mean(tPSTH(i,6,1200:2200),3)))); 
            else
                PsychEff(i) = ((AllConditions(i,4) / (AllConditions(i,4) + AllConditions(i,6))) - (AllConditions(i,5) / (AllConditions(i,5) + AllConditions(i,3)))); 
                NeuIdEff(i) = ((mean(tPSTH(i,4,1200:2200),3) / (mean(tPSTH(i,4,1200:2200),3) + mean(tPSTH(i,6,1200:2200),3))) - (mean(tPSTH(i,5,1200:2200),3) / (mean(tPSTH(i,5,1200:2200),3) + mean(tPSTH(i,3,1200:2200),3)))); 
            end
            ors (i)= GetProp(AllNeurons(i), FileType, StimulusType, 'or');
        end
end

figure(3250), clf, 
brucescatter(PsychEff, NeuIdEff, 1, [], filenamesforbruce); 
hold on, 
scatter(PsychEff(ors == 0 | ors == 180), NeuIdEff(ors == 0 | ors == 180),'r', 'filled');
scatter(PsychEff(abs(ors) ==90 | ors == 270), NeuIdEff(abs(ors) ==90 | ors == 270),'b', 'filled');

%%
rocs(rocs==-1)=0;
figure(918), scatter(rocs(:,1), rocs(:,2));
refline(0,0.5)

%% 

figure(19734), clf, hold on, 
brucescatter(rocs(:,1), rocs(:,2), 1, [], filenamesforbruce);
refline(0, 0.5);

%%
% figure, scatter(squeeze(mean(tPSTH(:,3,50:550),3)), squeeze(mean(tPSTH(:,5,50:550),3)), 'filled')
% hold on, scatter(squeeze(mean(tPSTH(:,4,50:550),3)), squeeze(mean(tPSTH(:,6,50:550),3)), 'filled')
% refline(1,0)
% refline(1,10)
% refline(1,-10)
% hold on, scatter(squeeze(mean(tPSTH(criteria,4,50:550),3)), squeeze(mean(tPSTH(criteria,6,50:550),3)))
% hold on, scatter(squeeze(mean(tPSTH(criteria,3,50:550),3)), squeeze(mean(tPSTH(criteria,5,50:550),3)))

%%
figure(19747), scatter(squeeze(mean(tPSTH(:,3,300:800),3)), squeeze(mean(tPSTH(:,5,300:800),3)), 'filled')
hold on,       scatter(squeeze(mean(tPSTH(:,4,300:800),3)), squeeze(mean(tPSTH(:,6,300:800),3)), 'filled')
refline(1,0)
refline(1,10)
refline(1,-10)
hold on, scatter(squeeze(mean(tPSTH(criteria,4,300:800),3)), squeeze(mean(tPSTH(criteria,6,300:800),3)))
hold on, scatter(squeeze(mean(tPSTH(criteria,3,300:800),3)), squeeze(mean(tPSTH(criteria,5,300:800),3)))

%%

figure(19834), clf, hold on, 
brucescatter(squeeze(mean(tPSTH(criteria,4,300:800),3)), squeeze(mean(tPSTH(criteria,6,300:800),3)), 1, [], filenamesforbruce(criteria))
brucescatter(squeeze(mean(tPSTH(criteria,3,300:800),3)), squeeze(mean(tPSTH(criteria,5,300:800),3)), 2, [], filenamesforbruce(criteria))

%%
figure(19844), clf, hold on, 
brucescatter(squeeze(mean(tPSTH(:,4,300:800),3)), squeeze(mean(tPSTH(:,6,300:800),3)), 1, [], filenamesforbruce)
brucescatter(squeeze(mean(tPSTH(:,3,300:800),3)), squeeze(mean(tPSTH(:,5,300:800),3)), 2, [], filenamesforbruce)


%%

figure(19934), clf, hold on, 
brucescatter(squeeze(mean(tPSTH(criteria,4,300:800),3))- squeeze(mean(tPSTH(criteria,6,300:800),3)), squeeze(mean(tPSTH(criteria,3,300:800),3))- squeeze(mean(tPSTH(criteria,5,300:800),3)), 2, [], filenamesforbruce(criteria))


%% relation between sd of gaussian fit to psychomtric function and the size of Id effect in BDID
c= 0;
for threshold = 0.2: -0.001:0  
    c = c + 1;
    criteria = ((abs(TI)>0.1) & (vs(2,:) >0) & (vs(2,:)<threshold)); 
    p = squeeze(mean(tPSTH(criteria,:,:)));
    IdE(c) = mean(p(7,550:end) - p(8,550:end));
    IdSE(c) = std(p(7,550:end) - p(8,550:end))/sqrt(size(p,2)-550);
    th(c) = threshold;
    nc(c) = sum(criteria);
end

figure(313), clf, hold on, 
scatter(th, IdE, 'filled');


figure(314), clf, hold on, 
scatter(nc, IdE, 'filled');

figure(315), clf, hold on,
errorbar(IdE, IdSE);

%%

% figure(17), scatter(smeb(criteria,1), smeb(criteria,7), 'filled'); title(''); xlabel('dx pref, bd pref'); ylabel('dx pref, bd null'); refline(1);
% figure(28), scatter(smeb(criteria,2), smeb(criteria,8), 'filled'); title(''); xlabel('dx null, bd pref'); ylabel('dx null, bd null'); refline(1);
% figure(36), scatter(smeb(criteria,3), smeb(criteria,6), 'filled'); title(''); xlabel('dx zero, bd pref'); ylabel('dx zero, bd null'); refline(1);
% 
% figure(110), clf, hold on
% scatter(smeb(criteria,1), smeb(criteria,7), 'r', 'filled');
% scatter(smeb(criteria,2), smeb(criteria,8), 'b', 'filled');
% scatter(smeb(criteria,3), smeb(criteria,6), 'g', 'filled');
% xlabel('bd pref'); ylabel('bd null'); 
% refline(1);
% title('');



figure(817), clf, brucescatter(smeb(criteria,1), smeb(criteria,7), 1, 7, filenamesforbruce(criteria)); title(''); xlabel('dx pref, bd pref'); ylabel('dx pref, bd null'); refline(1);
figure(828), clf, brucescatter(smeb(criteria,2), smeb(criteria,9), 3, 7, filenamesforbruce(criteria)); title(''); xlabel('dx null, bd pref'); ylabel('dx null, bd null'); refline(1);
figure(836), clf, brucescatter(smeb(criteria,3), smeb(criteria,6), 2, 7, filenamesforbruce(criteria)); title(''); xlabel('dx zero, bd pref'); ylabel('dx zero, bd null'); refline(1);



