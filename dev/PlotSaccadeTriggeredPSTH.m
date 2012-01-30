function [PSTHSacTrig, conditions, eb, es, varargout]= PlotSaccadeTriggeredPSTH(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, BinSize, Cumulative, PlotIt)
% Plots the Saccade Triggered PSTH for different conditions in the Experiment
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
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat');

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

try
    rdsPrefDir = PreferredRDSDirection(MonkeyName, NeuronNumber, ClusterName);
catch
    disp('rdsPrefDir FAILED HERE ! ! ');
    rdsPrefDir = Expt.Stimvals.or;
end

conditions = GetConditions(Expt, FileType, pD, rdsPrefDir);

SmoothingBinSize = 50;%50
SmthKernel = gausswin(SmoothingBinSize);
SmthKernel(SmoothingBinSize/2+1:end) = 0;

TrialLength = 500;
preSaccadeInterval = 100;
StartOffset = 500;

PSTHSacTrig = num2cell(zeros(length([Expt.Trials]),1));
PSTHSacTrigR= num2cell(zeros(length([Expt.Trials]),1));
shufTrials = randperm(length([Expt.Trials]));

scaling = Expt.Header.emtimes(end) / length(Expt.Header.emtimes) / 10;


%par
for tr = 1: length([Expt.Trials]),
    psth = [];
    saccade = [];
    scounter = 0;
    for sc = 1: length([Expt.Trials(tr).Saccades])
        if ((Expt.Trials(tr).Saccades(sc).start > StartOffset * 10) && (Expt.Trials(tr).Saccades(sc).start + TrialLength * 10 < Expt.Trials(tr).End))
          if (Expt.Trials(tr).Saccades(sc).size > 0.2)
            scounter = scounter + 1;
            for tt = 1:TrialLength,
                psth(tt) = sum([Expt.Trials(tr).Spikes]>=(tt*10 + Expt.Trials(tr).Saccades(sc).start - preSaccadeInterval * 10) & [Expt.Trials(tr).Spikes]<=(tt*10 + Expt.Trials(tr).Saccades(sc).start + BinSize*10 - preSaccadeInterval * 10)) * 1000 / BinSize;
            end
            psth = conv(psth, SmthKernel, 'same') ./ sum(SmthKernel);        
            PSTHSacTrig{tr, scounter} = psth;
            PSTHSacTrigR{shufTrials(tr), scounter} = psth;
            
            idata = Expt.Trials(tr).EyeData;
            saccade = idata(round(Expt.Trials(tr).Saccades(sc).start / scaling / 10): round(Expt.Trials(tr).Saccades(sc).end / scaling /10),:);
            %saccade = conv(saccade, SmthKernel, 'same') ./ sum(SmthKernel);        
            Saccades{tr, scounter} = saccade;
          end
        end
    end
end

eb = zeros(size(conditions,1), 0, TrialLength);
ebR = zeros(size(conditions,1), 0, TrialLength);
for cc = 1: size(conditions,1)
    for cndc = 1: size(conditions(cc,:),2),
        if (conditions(cc, cndc))
            for sc = 1 : size(PSTHSacTrig{cndc},1)
                eb(cc, size(eb,2) + 1, :) = mean(PSTHSacTrig{cndc}(sc, :),1);
                ebR(cc, size(ebR,2) + 1, :) = mean(PSTHSacTrigR{cndc}(sc, :),1);
            end
        end
    end
end


es = zeros(size(Saccades,1), 0, TrialLength);
AvgSaccade = zeros(size(Saccades,1), 4, 100);
s = [];
for iN = 1:size(Saccades,1)
    for sc = 1:size(Saccades(iN,:),2)
        if isempty(Saccades{iN,sc})
            break;
        else
            if (size(Saccades{iN,sc},1) > 5)
                s(sc,1:4,1:length(Saccades{iN,sc})) = Saccades{iN,sc}';
            end
        end
    end
    if (~isempty(s))
        if (size(s,1)>1)
            AvgSaccade(iN,1:4,1:size(s,3)) = squeeze(mean(s));
        else
            AvgSaccade(iN,1:4,1:size(s,3)) = squeeze(s);
        end
    end
end



for nanc = 1: size(eb,1)
    if sum(isnan(eb(nanc,:)))>1
        eb(nanc,:) = 0;
        se(nanc,:) = 0;
        disp('NAN detected!! here ############################## # # # # # ');
    end
end


warning off
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
        case 'DPI'
            figure(15121), h = plot(squeeze(mean(eb,2))');
            legend(h, GetLegends(FileType));
            refline(0);
        otherwise
            figure(19635),h = plot(eb');
            legend(h, GetLegends(FileType));
    end
    set(h, 'LineWidth', 2);

    set(gca, 'XGrid', 'on');
    xlim([0 500]);
    xtl = [-100, 0, 50, 100, 250, 500];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    title(strcat('NeuronID:',num2str(NeuronNumber), ' - BinSize:', num2str(BinSize), ' - PreferredCylinerRotationDir: ' , num2str(pD)));
    
    
end
warning on
