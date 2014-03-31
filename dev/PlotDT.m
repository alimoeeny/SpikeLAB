function [eb, se, values, pD]= PlotDT(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, PlotIt)
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

    pD = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
end

values = unique([Expt.Trials.(Expt.Stimvals.et)]);
[StartTime, FinishTime] = GetStartFinishTimes(FileType);

FRs = zeros(length([Expt.Trials]),1);
for tr = 1: size(FRs,1),
   FRs(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime) * 1000 / (FinishTime - StartTime);
end

for cc = 1: length(values)
    eb(cc,:) = mean(FRs([Expt.Trials(:).dx]==values(cc)));
    se(cc,:) =  std(FRs([Expt.Trials(:).dx]==values(cc))) / sqrt(sum([Expt.Trials(:).dx]==values(cc)));
end


for nanc = 1: size(eb,1)
    if sum(isnan(eb(nanc,:)))>1
        eb(nanc,:) = 0;
        se(nanc,:) = 0;
        disp('NAN detected!! here ############################## # # # # # ');
    end
end


if(PlotIt)
    figure, plot(eb), hold on, errorbar(eb,se);
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
%    xlim([-100 2500]);
%    xtl = [0, 200, 250, 700, 1200, 1700, 2200];
%    set(gca, 'XTick', xtl-BinSize/2);
%    set(gca, 'XTickLabel', {num2str((xtl-200)')});
    title(strcat('NeuronID:',num2str(NeuronNumber), ' - BinSize:', num2str(BinSize), ' - PreferredCylinerRotationDir: ' , num2str(pD)));    
end
