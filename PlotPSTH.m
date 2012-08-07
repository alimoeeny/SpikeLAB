function [PSTH, conditions, eb, varargout]= PlotPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, Cumulative, PlotIt, savetheplot)
% Plots the PSTH for different conditions in the Experiment
varargout{1} = [];
varargout{2} = [];

if(nargin==1)
    Expt = MonkeyName;
    FileType = 'NNN';
    BinSize = 50;
    PlotIt = 1;
    pD = PreferredCylinderRotationDirection(Expt);
    NeuronNumber = 0;
    Cumulative = 0;
else
    DataPath = GetDataPath();
    %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
%    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
 
    Neuron = load(filepath);
    Expt = Neuron.Expt;


    if strfind(FileType, 'RID')
         pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    %     if pD == -1 
    %         pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, 'DT', 0);
    %     end
    else
        pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
    end
end

try
    rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName);
catch
    disp('rdsPrefDir FAILED HERE ! ! ');
    rdsPrefDir = Expt.Stimvals.or;
end

conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);

SmoothingBinSize = 50;%50
SmthKernel = gausswin(SmoothingBinSize);
% make it causal
%SmthKernel(SmoothingBinSize/2+1:end) = 0;

StartTime = -2000; %1; %5500;
if(isfield(Expt.Trials, 'dur'))
    FinishTime = round(median([Expt.Trials(:).dur])) + 500;
else
    FinishTime = round(median([Expt.Trials(:).End] - [Expt.Trials(:).Start])) + 500;
end

PSTH = zeros(length([Expt.Trials]), round((FinishTime - StartTime + 1)/10));
%par
parfor tr = 1: size(PSTH,1),
    psth = [];
    for tt = 1:size(PSTH,2),
        psth(tt) = sum([Expt.Trials(tr).Spikes]>=(tt*10 + StartTime) & [Expt.Trials(tr).Spikes]<=(tt*10 + StartTime + BinSize*10)) * 1000 / BinSize;
    end
%     if Cumulative == 1
%         PSTH(tr,:) = cumsum(PSTH(tr,:));
%     end
    psth = conv(psth, SmthKernel, 'same') ./ sum(SmthKernel);        
    psths{tr} = psth;
end
for tr = 1: size(PSTH,1),
    PSTH(tr,:)= psths{tr};
end

for cc = 1: size(conditions,1)
    if sum(conditions(cc,:))>0
        eb(cc,:) = mean(PSTH(conditions(cc,:),:),1);
        se(cc,:) = std(PSTH(conditions(cc,:),:),[],1) / sqrt(sum(conditions(cc,:)));
    elseif cc > 1
        eb(cc,:) = 0;
        se(cc,:) = 0;
    end
end

if strcmpi(FileType, 'DID')
    startPointer = 750; %250;
    endPointer = 1500; %750;
    cuttoff = median(mean(PSTH(conditions(3,:), startPointer: endPointer),2));
    eDb(1,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) > cuttoff), :),1);
    eDb(2,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(1,:);
    eb(size(eb,1)+1,:) = eDb(2,:);
    conditions(size(conditions,1)+1,:)  = (conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff);

    cuttoff = median(mean(PSTH(conditions(4,:), startPointer: endPointer),2));
    eDb(3,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) > cuttoff), :),1);
    eDb(4,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(3,:);
    eb(size(eb,1)+1,:) = eDb(4,:);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff);

% ---------
    startPointer = 750;
    endPointer = 1500;
    cuttoff = median(mean(PSTH(conditions(3,:), startPointer: endPointer),2));
    eDb(5,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) > cuttoff), :),1);
    eDb(6,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(5,:);
    eb(size(eb,1)+1,:) = eDb(6,:);
    conditions(size(conditions,1)+1,:)  = (conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(3,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff);

    cuttoff = median(mean(PSTH(conditions(4,:), startPointer: endPointer),2));
    eDb(7,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) > cuttoff), :),1);
    eDb(8,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(7,:);
    eb(size(eb,1)+1,:) = eDb(8,:);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, startPointer: endPointer),2) <=cuttoff);

end


for nanc = 1: size(eb,1)
    if sum(isnan(eb(nanc,:)))>1
        eb(nanc,:) = 0;
        se(nanc,:) = 0;
        disp('NAN detected!! here ############################## # # # # # ');
    end
end




% choice probabilities
if (0) %////=============================
switch FileType
    case 'ABD'
        disp(ROCAUC(PSTH(conditions(1,:),1000:2000),PSTH(conditions(2,:),1000:2000))); 
    case {'DID', 'DIDB'}
        %roc1 = ROCAUC(PSTH(conditions(3,:),50:550),PSTH(conditions(5,:),50:550));
        %roc2 = ROCAUC(PSTH(conditions(6,:),50:550),PSTH(conditions(4,:),50:550));
        %disp(['First 500ms ROCs, Pref: ' num2str(roc1) ' , Null: ' num2str(roc2)]); 
        roc1 = ROCAUC(mean(PSTH(conditions(3,:),300:800),2),mean(PSTH(conditions(5,:),300:800),2));
        roc2 = ROCAUC(mean(PSTH(conditions(6,:),300:800),2),mean(PSTH(conditions(4,:),300:800),2));
        disp(['First 100-600ms ROCs, Pref: ' num2str(roc1) ' , Null: ' num2str(roc2)]); 
        a = zscore(mean(PSTH(conditions( 9,:)|conditions(10,:), 200:2100),2)); %1100:2100),2));
        b = zscore(mean(PSTH(conditions(11,:)|conditions(12,:), 200:2100),2)); %1100:2100),2));
        aa = [a(1 : sum(conditions(9,:))) ; b(1 : sum(conditions(11,:)))];
        bb = [a(sum(conditions(9,:))+1 : sum(conditions(9,:))+sum(conditions(10,:))) ; b(sum(conditions(11,:))+1 : sum(conditions(11,:))+sum(conditions(12,:)))];
        roc3 = ROCAUC(aa, bb);
        %main effect ROC 
        %roc4 = ROCAUC(PSTH(conditions(1,:),1200:2200),PSTH(conditions(2,:),1200:2200));
        [roc4, sig4, tmp] = ROCAUCSignificance(mean(PSTH(conditions(1,:),1200:2200),2),mean(PSTH(conditions(2,:),1200:2200),2));
        disp(['-------- --------- ---------- ------------ Main Effect sig:    ' , num2str(roc4), '   ===   ', num2str(sig4)]);
        roc5 = ROCAUC(PSTH(conditions(1,:),1000:2000),PSTH(conditions(2,:),1000:2000));
        %next to zero Stimulus AND choice
        roc6 = ROCAUC(PSTH(conditions(13,:),300:2000),PSTH(conditions(14,:),300:2000));
        % next to zero z scored
        a = zscore(mean(PSTH(conditions( 16,:)|conditions(17,:), 200:2100),2)); %1100:2100),2));
        b = zscore(mean(PSTH(conditions(18,:)|conditions(19,:), 200:2100),2)); %1100:2100),2));
        aa = [a(1 : sum(conditions(16,:))) ; b(1 : sum(conditions(18,:)))];
        bb = [a(sum(conditions(16,:))+1 : sum(conditions(16,:))+sum(conditions(17,:))) ; b(sum(conditions(18,:))+1 : sum(conditions(18,:))+sum(conditions(19,:)))];
        roc7 = ROCAUC(aa, bb);
        
        
        varargout{1}(1) = roc1;
        varargout{1}(2) = roc2;
        varargout{1}(3) = roc3;
        varargout{1}(4) = roc4;
        varargout{1}(5) = roc5;
        varargout{1}(6) = roc6;
        varargout{1}(7) = roc7;
    case 'BDID'
        %main effect ROC 
        roc4 = ROCAUC(PSTH(conditions(7,:),1200:2200),PSTH(conditions(8,:),1200:2200));
        varargout{1}(4) = roc4;
    case 'DPI' 
        eb(size(eb,1)+1, :) = eb(1,:) - eb(2,:);
        eb(size(eb,1)+1, :) = eb(3,:) - eb(4,:);
        eb(size(eb,1)+1, :) = eb(5,:) - eb(6,:);
end
end
%=========================================

warning off
if(PlotIt)
    switch upper(FileType) 
        case {'SRID', 'DRID'}
            figure,h = plot(eb(1:4,:)');
            legend('Pref Or Pref dx','Pref Or null dx','Null Or Pref dx ','null Or null dx');
        case 'DID'
            disp(['Max dx:', num2str(max([Expt.Trials(:).dx])), ' - Ids are: ' , num2str(unique([Expt.Trials(:).Id]))]);
            trialcombos = [1,2];% [1, 2, 15]; %[3,4,5,6]; %[1:size(eb,1)]; % [1,2,7,8];
            hp = figure(18719); h = plot(eb(trialcombos,:)'); % plot(eb'); %plot(eb([2,4,6],:)');
            leg = GetLegends(FileType);
            legend(h, leg(trialcombos));
        case 'TWO'
            trialcombos = [3,6]; %[1:size(eb,1)];
            hp = figure(13719); clf; h = plot(eb(trialcombos,:)');
            leg = GetLegends(FileType);
            legend(h, leg(trialcombos));
        case 'BDID'
            trialcombos = [7,8];
            hp = figure(17819); h = plot(eb(trialcombos,:)');
            leg = GetLegends(FileType);
            legend(h, leg(trialcombos));
        case 'DPI'
            figure(16919), h = plot(eb');
            legend(h, GetLegends(FileType));
            refline(0);
        case 'DTRW'
            disp(['Max dx:', num2str(max([Expt.Trials(:).dx])), ' - Rewards are: ' , num2str(unique([Expt.Trials(:).rw]))]);
            trialcombos = [1:size(eb,1)];% [1, 2, 15]; %[3,4,5,6]; %[1:size(eb,1)]; % [1,2,7,8];
            hp = figure(12719); h = plot(eb(trialcombos,:)'); % plot(eb'); %plot(eb([2,4,6],:)');
            leg = GetLegends(FileType);
            legend(h, leg(trialcombos));
        otherwise
            figure(19635),h = plot(eb');
            legend(h, GetLegends(FileType));
    end
    set(h, 'LineWidth', 2);

    set(gca, 'XGrid', 'on');
    xlim([-100 2500]);
    xtl = [0, 200, 250, 700, 1200, 1700, 2200];
    set(gca, 'XTick', xtl-BinSize/2);
    set(gca, 'XTickLabel', {num2str((xtl-200)')});
    title(strcat('NeuronID:',num2str(NeuronNumber), ' - BinSize:', num2str(BinSize), ' - PreferredCylinerRotationDir: ' , num2str(pD)));
    
    if(savetheplot)
      print(hp, '-dpsc', '-r300', '-zbuffer', ['../figs/PSTH', '-', date, '-', filename, '.eps']);
      
      figure(hp+110); clf;
      px = PlotExpt(filepath);
      print(px.fig, '-dpsc', '-r300', '-zbuffer', ['../figs/PlotExpt', '-', date, '-', filename, '.eps']);
      
    end
    
%     hold on,
% %    errorbar(eb', se');
%     xstep = 100;
%     xfake = 1:xstep:size(eb,2);
%     for xfk = 2:size(eb,1)
%         xfake = [xfake;1:xstep:size(eb,2)] ;
%     end
%     if strfind(FileType, 'RID')
%        errorbar(xfake(:,1:4), eb(1:4,[1:xstep:end])', se(1:4,[1:xstep:end])', 'x', 'LineWidth',2)   
%     else
%        errorbar(xfake', eb(:,[1:xstep:end])', se(:,[1:xstep:end])', 'x', 'LineWidth',2)   
%     end
    
end
warning on
