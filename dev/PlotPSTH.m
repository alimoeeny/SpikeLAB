function [PSTH, conditions, eb, varargout]= PlotPSTHIdisp(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, Cumulative, PlotIt)
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
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');

    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
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

conditions = GetConditions(Expt, FileType, pD);

StartTime = -2000; %1; %5500;
if(isfield(Expt.Trials, 'dur'))
    FinishTime = round(median([Expt.Trials(:).dur])) + 500;
else
    FinishTime = round(median([Expt.Trials(:).End] - [Expt.Trials(:).Start])) + 500;
end

PSTH = zeros(length([Expt.Trials]), round((FinishTime - StartTime + 1)/10));
parfor tr = 1: size(PSTH,1),
    psth = [];
    for tt = 1:size(PSTH,2),
        psth(tt) = sum([Expt.Trials(tr).Spikes]>=(tt*10 + StartTime) & [Expt.Trials(tr).Spikes]<=(tt*10 + StartTime + BinSize*10)) * 1000 / BinSize;
    end
%     if Cumulative == 1
%         PSTH(tr,:) = cumsum(PSTH(tr,:));
%     end
    psths{tr} = psth;
end
for tr = 1: size(PSTH,1),
    PSTH(tr,:)= psths{tr};
end

for cc = 1: size(conditions,1)
    if sum(conditions(cc,:))>0
        eb(cc,:) = mean(PSTH(conditions(cc,:),:),1);
        se(cc,:) = std(PSTH(conditions(cc,:),:),[],1) / sqrt(sum(conditions(cc,:)));
    else
        eb(cc,:) = 0;
        se(cc,:) = 0;
    end
end

if strcmpi(FileType, 'DID')
    cuttoff = median(mean(PSTH(conditions(3,:), 250: 750),2));
    eDb(1,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, 250: 750),2) > cuttoff), :),1);
    eDb(2,:)      = mean(PSTH((conditions(3,:)' & mean(PSTH(:, 250: 750),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(1,:);
    eb(size(eb,1)+1,:) = eDb(2,:);
    conditions(size(conditions,1)+1,:)  = (conditions(3,:)' & mean(PSTH(:, 250: 750),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(3,:)' & mean(PSTH(:, 250: 750),2) <=cuttoff);

    cuttoff = median(mean(PSTH(conditions(4,:), 250: 750),2));
    eDb(3,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, 250: 750),2) > cuttoff), :),1);
    eDb(4,:)      = mean(PSTH((conditions(4,:)' & mean(PSTH(:, 250: 750),2) <=cuttoff), :),1);
    varargout{2} = eDb;
    eb(size(eb,1)+1,:) = eDb(1,:);
    eb(size(eb,1)+1,:) = eDb(2,:);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, 250: 750),2) >cuttoff);
    conditions(size(conditions,1)+1,:) = (conditions(4,:)' & mean(PSTH(:, 250: 750),2) <=cuttoff);
end


for nanc = 1: size(eb,1)
    if sum(isnan(eb(nanc,:)))>1
        eb(nanc,:) = 0;
        se(nanc,:) = 0;
        disp('NAN detected!! here ############################## # # # # # ');
    end
end

% choice probabilities
switch FileType
    case 'ABD'
        disp(ROCAUC(PSTH(conditions(1,:),1000:2000),PSTH(conditions(2,:),1000:2000))); 
    case {'DID', 'DIDB'}
        %roc1 = ROCAUC(PSTH(conditions(3,:),50:550),PSTH(conditions(5,:),50:550));
        %roc2 = ROCAUC(PSTH(conditions(6,:),50:550),PSTH(conditions(4,:),50:550));
        %disp(['First 500ms ROCs, Pref: ' num2str(roc1) ' , Null: ' num2str(roc2)]); 
        roc1 = ROCAUC(PSTH(conditions(3,:),300:600),PSTH(conditions(5,:),300:800));
        roc2 = ROCAUC(PSTH(conditions(6,:),300:600),PSTH(conditions(4,:),300:800));
        disp(['First 100-600ms ROCs, Pref: ' num2str(roc1) ' , Null: ' num2str(roc2)]); 
        a = zscore(mean(PSTH(conditions( 9,:)|conditions(10,:),1100:2100),2));
        b = zscore(mean(PSTH(conditions(11,:)|conditions(12,:),1100:2100),2));
        aa = [a(1 : sum(conditions(9,:))) ; b(1 : sum(conditions(11,:)))];
        bb = [a(sum(conditions(9,:))+1 : sum(conditions(9,:))+sum(conditions(10,:))) ; b(sum(conditions(11,:))+1 : sum(conditions(11,:))+sum(conditions(12,:)))];
        roc3 = ROCAUC(aa, bb);
        varargout{1}(1) = roc1;
        varargout{1}(2) = roc2;
        varargout{1}(3) = roc3;
end

if(PlotIt)
    switch upper(FileType) 
        case {'SRID', 'DRID'}
            figure,h = plot(eb(1:4,:)');
            legend('Pref Or Pref dx','Pref Or null dx','Null Or Pref dx ','null Or null dx');
        case 'DID'
            disp(['Max dx:', num2str(max([Expt.Trials(:).dx])), ' - Ids are: ' , num2str(unique([Expt.Trials(:).Id]))]);
            figure(18719),h = plot(eb'); %plot(eb([2,4,6],:)');
            legend(h, GetLegends(FileType));
            % %% figure,h = plot(eb([1,2,7,8],:)');
            % %% legend(h, 'Preferred Id', 'Null Id', 'Pref dx', 'Null dx');
        case 'TWO'
            figure,h = plot(eb');
            legend('Pref bd Pref dx', 'Pref bd null dx','Pref bd ZERO dx', 'Pref bd ZERO dx and correct response', 'null bd ZERO dx', 'null bd ZERO dx correct response');
        case 'BDID'
            figure(17819),h = plot(eb(7:8,:)');
            legend(h, GetLegends(FileType));
            %legend('Pref bd Pref Id', 'Pref bd null Id', 'Pref bd pref Id and correct response', 'Pref bd Pref Id wrong response', 'null bd null Id correct response', 'null bd null Id wrong response', 'Preff Id', 'Null Id');
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