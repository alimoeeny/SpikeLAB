clear;
%clc;

ShowSingleCellSDFs = 0; % 0 or 1
[AllNeurons, FileType, StimulusType] = loadAllNeurons4('DID');
%Prep
DataPath = GetDataPath();
BinSize = 1;%50;
SmoothingBinSize = 1;%50
SmthKernel = gausswin(SmoothingBinSize);


filenamesforbruce = {};
TI=[];

% AllNeurons =  SelectByMonkey(AllNeurons, 'ic');
% AllNeurons =  SelectByMonkey(AllNeurons, 'dae');
% disp('  O N E   M O N K E Y   A T  A  T I M E ');


%AllNeurons = AllNeurons([3,4,6,7, 11, 13, 15, 20, 21, 22,23]);

SacBars = {};
%par
for iN= 1: length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));

    p = []; c = []; eb = []; jnk1 = []; jnk2 = []; xC =[]; c1 =[]; ebX =[]; 
    pSacTrig=[]; cSacTrig=[]; ebSacTrig=[]; ebSacTrigMean = [];
    
    switch FileType 
        case {'ABD', 'DID', 'DRID', 'SRID', 'DIDB'} 
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, [], 1);
            %if (abs(TI(iN))>0.1)
                %[xC ,c1, ebX] = PlotSTAutoCorr(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, ShowSingleCellSDFs);
                [p, c, eb, jnk1, jnk2] = PlotPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
                %[c, eb] = PlotISISdf(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, 1);
                slp{iN} = PsychSlop(AllNeurons(iN), StimulusType, FileType, 'dx'); 
                %disp(['Ti; ', num2str(TI(iN)), ' - Slop: ' , num2str(slp(iN).fit(2))]);
            %end
        case 'BDID'
            [p, c, eb, jnk1, jnk2] = PlotPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
            slp{iN} = PsychSlop(AllNeurons(iN), StimulusType, FileType, 'bd'); 
            %disp(['Ti: ', num2str(TI(iN)), ' - Slop: ' , num2str(slp(iN).fit(2))]);
        case 'TWO'
            [p, c, eb] = PlotPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
        case 'DPI'
            ti      = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, [], 1);
            %TIcyldx(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, 'DT',     [], 1);
            %[p, c, eb] = PlotPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            [pSacTrig, cSacTrig, ebSacTrig] = PlotSaccadeTriggeredPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, 0, ShowSingleCellSDFs);
            %[pSacM, pSacSe, pSacPrj] = PlotSaccades(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType);
            %SacBars{iN} = [pSacM , pSacSe, pSacPrj];
            
            %PIS(iN,:) = PursuitIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, ti);
            %dxEffectPeak(iN) = dxEffectMax(MonkeyName, NeuronNumber, ClusterName, 'rds');
            TI(iN) = ti;
    end
    
    pD = -10;
    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    orS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    
    if (strcmp(FileType, 'DRID') || strcmp(FileType, 'SRID') || ((strcmpi(FileType, 'DPI')) && (strcmpi(StimulusType, 'rds'))))
        cellValid  = 1; vScore = 1;
    else
        [cellValid , vScore] = CellValidity([], FileType, pD, MonkeyName, NeuronNumber, ClusterName, FileType);
    end
    
    
    for iST = 1: size(ebSacTrig,1),
        ebSacTrigMean(iST,:) = squeeze(mean(ebSacTrig(iST,squeeze(mean(ebSacTrig(iST,:,:),3))'>0,:)));
    end
          
    ebSacTrig = []; pSacTrig = []; cSacTrig = []; % We dont need these for now
    PSTHs{iN} = {p , c, eb, cellValid, vScore, pD, jnk1, jnk2, xC, ebX, pSacTrig, cSacTrig, ebSacTrig, ebSacTrigMean};
    
    filenamesforbruce{iN} = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    
    binoc{iN} = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 've'); %Expt.Stimvals.ve;
    %save(['~/Desktop/matlab.', num2str(iN),'.mat']);
end


%%   
save ~/Desktop/matlab.mat -v7.3

%%
validCells = []; for i = 1: length(PSTHs), if~isempty(PSTHs{i}), validCells(i) = [PSTHs{i}{4}]; end; end
AllConditions = [];

for i = 1 : length(PSTHs) % length(AllNeurons),
    if (~isempty(SacBars))
        if (~isempty(SacBars{i}))
            SacMeanVects(i,:,:,:,:) = SacBars{i};
            SacSEMVects (i,:,:,:,:) = SacBars{i};
        end
    end
    if (~isempty(PSTHs{i}) && ~isempty(PSTHs{i}{14}))
        tPSTHSacTrigMean(i, :, :) = PSTHs{i}{14};
    end
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

        
        tPSTHSacTrigMean(i, :, :) = PSTHs{i}{14};
        
        cntr = 0;
        for ebi = 1:size(eb,1)
            cntr = cntr + 1;
            if(sum(eb(ebi,:))>0)
                tPSTH(i, cntr,1:size(eb(cntr,:),2)) = eb(ebi,:);        
                if(ebi<=size(ebX,1))
                    txCor(i, cntr,1:size(ebX(cntr,:),2))= ebX(ebi,:);
                end
                
                nidp(i,1) = (sum(PSTHs{i}{2}(3,:))); 
                nidp(i,2) = (sum(PSTHs{i}{2}(4,:)));
                nidp(i,3) = (sum(PSTHs{i}{2}(5,:))); 
                nidp(i,4) = (sum(PSTHs{i}{2}(6,:)));
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
    case 'DPI'
        criteria = (TI>0.1 | TI<-0.1);% & validCells;
end


%criteria = criteria & ([binoc{:}]>5.07);
%criteria = criteria & (TI<0);
%figure, plot(cumsum(criteria));

%% 

% for n = 1:length(PIS)
%   disp(orS(n));
%   if ((orS(n)>0) & (orS(n)<=180)) 
%     if TI(n)>0 
%         % we are good
%     else
%         tmp = PIS(n, 2); PIS(n, 2) = PIS(n,3); PIS(n, 3) = tmp;
%         tmp = PIS(n, 5); PIS(n, 5) = PIS(n,6); PIS(n, 6) = tmp;
%         tmp = PIS(n, 8); PIS(n, 8) = PIS(n,9); PIS(n, 9) = tmp;
%     end
%   else
%     PIS(n,:) = - PIS(n,:);
%     if TI(n)>0 
%         % we are good
%     else
%         tmp = PIS(n, 2); PIS(n, 2) = PIS(n,3); PIS(n, 3) = tmp;
%         tmp = PIS(n, 5); PIS(n, 5) = PIS(n,6); PIS(n, 6) = tmp;
%         tmp = PIS(n, 8); PIS(n, 8) = PIS(n,9); PIS(n, 9) = tmp;
%     end
%   end
% end

% JUNK figure(1919), scatter(PIS(:,1) ./ PIS(:,7), PIS(:,4) ./ PIS(:,7), 'filled'); refline(1), refline(0,0), reflinexy(10,1), reflinexy(1,10), reflinexy(0,10);

figure(1919), scatter(PIS(:,1) ./ PIS(:,7), PIS(:,4) ./ PIS(:,7), 'filled'); refline(1), refline(0,0), reflinexy(10,1), reflinexy(1,10), reflinexy(0,10);
figure(1918), scatter(PIS(:,11) ./ PIS(:,13), PIS(:,12) ./ PIS(:,13), 'filled'); refline(1), refline(0,0), reflinexy(10,1), reflinexy(1,10), reflinexy(0,10);
figure(1215), clickscatter(PIS(:,15), dxEffectPeak, 6, 6, filenamesforbruce ); refline(1)

%%
figure(1920), clf, hold on,
scatter(PIS(:,11), PIS(:,12), 'filled'); 
scatter(PIS(:,11), PIS(:,13), 'filled'); 
refline(1), refline(0,0), reflinexy(10,1), reflinexy(1,10), reflinexy(0,10);

%%
if strcmp(FileType,'TWO')
    PopPSTH = mean(tPSTH(criteria,:,:));
else
    PopPSTH = squeeze(mean(tPSTH(criteria,:,:)));
    PopPSTHse = std(tPSTH(criteria,:,:))/sqrt(size(tPSTH,1));
    if exist('txCor')
        PopxCor = squeeze(mean(txCor(criteria,:,:)));
    end
end

%%
if exist('tPSTHSacTrigMean')
  PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean(criteria,:,:)));
end


%%

w_mode = 'trial_count';%'trial_count'; %'pref_init_trial';%'pref_trial'; %'pref'; %'trial_count';
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
            switch FileType
                case 'DID'
                    prefMean = mean(squeeze(tPSTH(wi, 1, :)));
                case 'BDID'
                    prefMean = mean(squeeze(tPSTH(wi, 7, :)));
            end
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = squeeze(tPSTH(wi,wc,:))/prefMean;
            end
        end
        WeightedPopPSTH = squeeze(mean(WtPSTH(criteria,:,:)));
    case 'pref_init'
        for wi =1: size(tPSTH,1)
            switch FileType
                case 'DID'
                    prefMean = mean(squeeze(tPSTH(wi, 1, 250:750)));
                case 'BDID'
                    prefMean = mean(squeeze(tPSTH(wi, 7, 250:750)));
            end
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = squeeze(tPSTH(wi,wc,:))/prefMean;
            end
        end
        WeightedPopPSTH = squeeze(mean(WtPSTH(criteria,:,:)));    
    case 'pref_trial'
        for wi =1: size(tPSTH,1)
            switch FileType
                case 'DID'
                    prefMean = mean(squeeze(tPSTH(wi, 1, :)));
                case 'BDID'
                    prefMean = mean(squeeze(tPSTH(wi, 7, :)));
            end
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
    case 'pref_init_trial'
        for wi =1: size(tPSTH,1)
            switch FileType
                case 'DID'
                    prefMean = mean(squeeze(tPSTH(wi, 1, 250:750)));
                case 'BDID'
                    prefMean = mean(squeeze(tPSTH(wi, 7, 250:750)));
            end
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
    case 'norm_full'
        for wi =1: size(tPSTH,1)
            switch FileType
                case 'DID'
                    prefMean = mean(squeeze(tPSTH(wi, 15, :)));
                case 'BDID'
                    prefMean = mean(squeeze(tPSTH(wi, 15, :)));
            end
            for wc = 1: size(AllConditions,2)
                WtPSTH(wi,wc,:) = squeeze(tPSTH(wi,wc,:))/prefMean;
            end
        end
        WeightedPopPSTH = squeeze(mean(WtPSTH(criteria,:,:)));

end



%%
switch FileType
    case 'SRID'
        figure(142);
    case 'DRID'
        figure(143);
    case 'DID'
        figure(162);
    case 'BDID'
        figure(136);  
    case 'DPI'
        figure(188);
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

%%
switch FileType
    case 'SRID'
        figure(642);
    case 'DRID'
        figure(643);
    case 'DID'
        figure(662);
    case 'BDID'
        figure(636);  
    case 'DPI'
        figure(689);
    otherwise
        figure(625);
end
clf, hold on,  
h = plot(squeeze(PopPSTHSacTrigMean)');

set(h, 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);

%errorbar(squeeze(PopPSTH)', squeeze(PopPSTHse)');

sum(criteria)



%% 

figure(2502),
clf, hold on,

h = plot(squeeze(mean(tPSTHSacTrigMean(PIS(:,10)<0,[41:42 71:76],:)))');
set(h, 'LineWidth', 2);
h = plot(squeeze(mean(tPSTHSacTrigMean(PIS(:,10)>0,[41:42 71:76],:)))');
set(h, 'LineWidth', 4);
h = plot(squeeze(mean(tPSTHSacTrigMean(:,[41:42 71:76],:)))');
set(h, 'LineWidth', 6);

set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);


%%
%load ~/Desktop/matlab.cylDPI.new.mat
figure(8392), clf , hold on,
a(:,1) = squeeze(mean(tPSTHSacTrigMean(TI>0,41,:)));
a(:,2) = squeeze(mean(tPSTHSacTrigMean(TI<0,41,:)));
a(:,3) = squeeze(mean(tPSTHSacTrigMean(TI>0,42,:)));
a(:,4) = squeeze(mean(tPSTHSacTrigMean(TI<0,42,:)));
h = plot(mean(a(:,[1,4]),2));
h = plot(mean(a(:,[2,3]),2));

%%
for i = 1: length(AllNeurons),
    PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean(1:i,41:end,:)));
    figure(007007), clf, hold on,
    h = plot(squeeze(PopPSTHSacTrigMean)');
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
    xlim([0 500]);
    xtl = [-100, 0, 50, 100, 250, 500];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    legend(h, GetLegends(FileType));
    title(FileType);
end



%%
load ~/Desktop/matlab.cylDPI.new.mat
figure(8392), clf , hold on,
a(:,1) = squeeze(mean(tPSTHSacTrigMean(TI>0,41,:)));
a(:,2) = squeeze(mean(tPSTHSacTrigMean(TI<0,41,:)));

a(:,3) = squeeze(mean(tPSTHSacTrigMean(TI>0,42,:)));
a(:,4) = squeeze(mean(tPSTHSacTrigMean(TI<0,42,:)));

h = plot(mean(a(:,[1]),2));
h = plot(mean(a(:,[2]),2));
h = plot(mean(a(:,[3]),2));
h = plot(mean(a(:,[4]),2));
h = plot(mean(a(:,[1,4]),2));
h = plot(mean(a(:,[2,3]),2));

set(h, 'LineWidth',2);
set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);

%%

mtn = 5; % minimum number of trials required
ng = []; SDFa1 = []; SDFa2 = []; SDFa3 = []; SDFb1 = []; SDFb2 = []; SDFb3 = []; SDF4_1 = []; SDF4_2 = [];
for i = 1:length(PSTHs)
    cnd1 = ((nidp(i,1)>=mtn) && (nidp(i,3)>=mtn));
    cnd2 = ((nidp(i,2)>=mtn) && (nidp(i,4)>=mtn));
    if ~(cnd1 || cnd2), 
        disp([num2str(i) AllNeurons(i)]);
        ng(i) = 1;
        SDF4_1 = [SDF4_1; squeeze(tPSTH(i,1,:))']; % Id pref
        SDF4_2 = [SDF4_2; squeeze(tPSTH(i,2,:))']; % Id null
    else if (cnd1 && cnd2), 
            ng(i) = 4;
            SDFa1 = [SDFa1; squeeze(tPSTH(i,3,:))']; %correct pref
            SDFa2 = [SDFa2; squeeze(tPSTH(i,4,:))']; %correct null
            SDFa3 = [SDFa3; squeeze(tPSTH(i,5,:))']; %wrong pref

            SDFb1 = [SDFb1; squeeze(tPSTH(i,3,:))']; %correct pref
            SDFb2 = [SDFb2; squeeze(tPSTH(i,4,:))']; %correct null
            SDFb3 = [SDFb3; squeeze(tPSTH(i,6,:))']; %wrong null

        else if (cnd1 && ~cnd2)
            ng(i) = 2;
            SDFa1 = [SDFa1; squeeze(tPSTH(i,3,:))']; %correct pref
            SDFa2 = [SDFa2; squeeze(tPSTH(i,4,:))']; %correct null
            SDFa3 = [SDFa3; squeeze(tPSTH(i,5,:))']; %wrong pref
   
            else
                ng(i) = 3;
                SDFb1 = [SDFb1; squeeze(tPSTH(i,3,:))']; %correct pref
                SDFb2 = [SDFb2; squeeze(tPSTH(i,4,:))']; %correct null
                SDFb3 = [SDFb3; squeeze(tPSTH(i,6,:))']; %wrong null
                
            end
        end
        
    end
end


figure(876), hist(ng);
%%

figure(888); clf, hold on,
plot(mean([SDFa1; SDFb1])', 'r', 'LineWidth', 2);
plot(mean([SDFa2; SDFb2])', 'b', 'LineWidth', 2);
plot(mean(SDFa3)', 'm', 'LineWidth', 2);
plot(mean(SDFb3)', 'c', 'LineWidth', 2);
plot(mean(SDF4_1)', 'k');
plot(mean(SDF4_2)', 'k');

set(gca, 'XGrid', 'on');
xlim([100 2200]);
xtl = [0, 50, 500, 1000, 1500, 2000];
set(gca, 'XTick', xtl+200-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);


%% Sliding Delta


%colors= 'rgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmyk';
switch FileType
    case {'BDID'}
        figure(7231), clf
        hold on,
        [vss, idx] = sort(vs);
        for i = 1:length(vs)-8
            %criteria = ((abs(TI)>0.1) & (vs >thresholds(i-1)) & (vs<thresholds(i)));
            criteria = idx(i:i+8);
            psq = squeeze(mean(tPSTH(criteria,:,:)))';
            if (abs(TI(i))>0.1)
                h = plot((psq(:,7) - psq(:,8)) ./ (psq(:,7) + psq(:,8)), 'Color', [(length(vs)-i)/400 (length(vs)-i)/150 (length(vs)-i)/200]);
            end
            d1(i) = mean(psq(200:800,7) - psq(200:800,8));
            d2(i) = mean(psq(800:2200,7) - psq(800:2200,8));
    
            set(h, 'LineWidth', 2);
        end
    case {'DID'}
        figure(7221), clf
        hold on,
        [vss, idx] = sort(vs);
        for i = 1:length(vs)-8
            %criteria = ((abs(TI)>0.1) & (vs >thresholds(i-1)) & (vs<thresholds(i)));
            criteria = idx(i:i+8);
            PopPSTH = squeeze(mean(tPSTH(criteria,:,:)));
            psq = squeeze(PopPSTH)';
            h = plot((psq(:,1) - psq(:,2)) ./ (psq(:,1) + psq(:,2)), 'Color', [(length(vs)-i)/100 (length(vs)-i)/200 (length(vs)-i)/400]);
            d1(i) = mean(psq(200:800,1) - psq(200:800,2));
            d2(i) = mean(psq(800:2200,1) - psq(800:2200,2));
    
            set(h, 'LineWidth', 2);
        end
end
        
set(gca, 'XGrid', 'on');
xlim([100 2200]);
xtl = [0, 50, 500, 1000, 1500, 2000];
set(gca, 'XTick', xtl+200-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});

refline(0,0)

figure(7222), plot(d1, 'r')
hold on, plot(d2, 'b')


%% Delta
colors= 'rgbcmyk';
switch FileType
    case {'BDIDB'}
        figure(723), clf
        hold on,
        thresholds = [0,0.01, 0.015, 0.1];

        for i = 2:length(thresholds)
            criteria = ((abs(TI)>0.1) & (vs >thresholds(i-1)) & (vs<thresholds(i)));
            PopPSTH = squeeze(mean(tPSTH(criteria,:,:)));
            psq = squeeze(PopPSTH)';
            h = plot(psq(:,7) - psq(:,8), colors(i-1));
            disp(sum(criteria))
            set(h, 'LineWidth', 2);
        end
    case {'DID'}
        figure(722), clf
        hold on,
        thresholds = [0,0.01, 0.02, 0.1];

        for i = 2:length(thresholds)
            criteria = ((abs(TI)>0.1) & (vs >thresholds(i-1)) & (vs<thresholds(i)));
            PopPSTH = squeeze(mean(tPSTH(criteria,:,:)));
            psq = squeeze(PopPSTH)';
            h = plot(psq(:,1) - psq(:,2), colors(i-1));
            disp(sum(criteria))
            set(h, 'LineWidth', 2);
        end
end
        
set(gca, 'XGrid', 'on');
xlim([100 2200]);
xtl = [0, 50, 500, 1000, 1500, 2000];
set(gca, 'XTick', xtl+200-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});

refline(0,0)

%% 

figure(456), clf,

clickscatter(rocs(:,4), rocs(:,3), 1, [], filenamesforbruce)


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
clickscatter(PsychEff, NeuIdEff, 1, [], filenamesforbruce); 
hold on, 
scatter(PsychEff(ors == 0 | ors == 180), NeuIdEff(ors == 0 | ors == 180),'r', 'filled');
scatter(PsychEff(abs(ors) ==90 | ors == 270), NeuIdEff(abs(ors) ==90 | ors == 270),'b', 'filled');

%%
rocs(rocs==-1)=0;
figure(918), scatter(rocs(:,1), rocs(:,2));
refline(0,0.5)

%% 

figure(19734), clf, hold on, 
clickscatter(rocs(:,1), rocs(:,2), 1, [], filenamesforbruce);
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
clickscatter(squeeze(mean(tPSTH(criteria,4,300:800),3)), squeeze(mean(tPSTH(criteria,6,300:800),3)), 1, [], filenamesforbruce(criteria))
clickscatter(squeeze(mean(tPSTH(criteria,3,300:800),3)), squeeze(mean(tPSTH(criteria,5,300:800),3)), 2, [], filenamesforbruce(criteria))

%%
figure(19844), clf, hold on, 
clickscatter(squeeze(mean(tPSTH(:,4,300:800),3)), squeeze(mean(tPSTH(:,6,300:800),3)), 1, [], filenamesforbruce)
clickscatter(squeeze(mean(tPSTH(:,3,300:800),3)), squeeze(mean(tPSTH(:,5,300:800),3)), 2, [], filenamesforbruce)


%%

figure(19934), clf, hold on, 
clickscatter(squeeze(mean(tPSTH(criteria,4,300:800),3))- squeeze(mean(tPSTH(criteria,6,300:800),3)), squeeze(mean(tPSTH(criteria,3,300:800),3))- squeeze(mean(tPSTH(criteria,5,300:800),3)), 2, [], filenamesforbruce(criteria))


%% relation between sd of gaussian fit to psychomtric function and the size of Id effect in BDID
c= 0;
for threshold = 0.2: -0.001:0  
    c = c + 1;
    criteria = ((abs(TI)>0.1) & (vs >0) & (vs<threshold)); 
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



figure(817), clf, clickscatter(smeb(criteria,1), smeb(criteria,7), 1, 7, filenamesforbruce(criteria)); title(''); xlabel('dx pref, bd pref'); ylabel('dx pref, bd null'); refline(1);
figure(828), clf, clickscatter(smeb(criteria,2), smeb(criteria,9), 3, 7, filenamesforbruce(criteria)); title(''); xlabel('dx null, bd pref'); ylabel('dx null, bd null'); refline(1);
figure(836), clf, clickscatter(smeb(criteria,3), smeb(criteria,6), 2, 7, filenamesforbruce(criteria)); title(''); xlabel('dx zero, bd pref'); ylabel('dx zero, bd null'); refline(1);


%% DPI index

for iN = 1: length(AllNeurons)
    [a,x, y, z] = ttest(PSTHs{1,iN}{3}(1,:) ,PSTHs{1,iN}{3}(2,:), 0.001);
    b = mean(PSTHs{1,iN}{3}(1,:) - PSTHs{1,iN}{3}(2,:)) / mean(PSTHs{1,iN}{3}(2,:));
    p = anova1([PSTHs{1,iN}{3}(1,:); PSTHs{1,iN}{3}(2,:); PSTHs{1,iN}{3}(3,:)]', [], 'off');
    
    conditions = PSTHs{1,iN}{2};
    
    d =  mean([mean(PSTHs{1,iN}{1}(conditions(3,:),:),2) ; mean(PSTHs{1,iN}{1}(conditions(4,:),:),2)]) - ...
         mean([mean(PSTHs{1,iN}{1}(conditions(5,:),:),2) ; mean(PSTHs{1,iN}{1}(conditions(6,:),:),2)]) ;
    e =  ( mean(mean(PSTHs{1,iN}{1}(conditions(1,:),:),2)) - mean(mean(PSTHs{1,iN}{1}(conditions(2,:),:),2)) );
    
    %d = mean( (PSTHs{1,iN}{3}(3,:) + PSTHs{1,iN}{3}(4,:)) - (PSTHs{1,iN}{3}(5,:) + PSTHs{1,iN}{3}(6,:)) ) * 0.5 ;
    %r = mean(  PSTHs{1,iN}{3}(1,:) - PSTHs{1,iN}{3}(2,:)) / d;
    
    disp([abs(TI(iN))>0.1,a, b, p, e, e / d]);
    aa = a;
    pp(iN) = p;
    bb(iN) = b;
    ee(iN) = e;
    dd(iN) = d;
end

disp(mean(ee) / mean(dd));


%%


figure(1312), 
clf, hold on,
clickscatter(SacMeanVects(:,41,5), SacMeanVects(:,42,5), 6, 6, filenamesforbruce);
clickscatter(SacMeanVects(:,71,5), SacMeanVects(:,42,5), 7, 7, filenamesforbruce);
   






