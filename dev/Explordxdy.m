% dxdy

clear, clc

cd /Users/moeenya/Dropbox/Projects/SpikeLAB/dev

DataPath = GetDataPath();
doitsquare = 1;

load ../AllDxDyNeurons.mat
AllNeurons = AllDxDyNeurons;
clear AllDxDyNeurons
FileType = 'DXDY';
StimulusType = 'rds';
StartTime  = 500;
FinishTime = 5000;

gfs=1024;
gfilter = zeros(gfs, gfs);
gw = 8;
gfilter(gfs/2 - gw/2 + 1:gfs/2+gw/2,gfs/2 - gw/2+1:gfs/2+gw/2) = gausswin(gw)* gausswin(gw)';


%%
for iN  = 1:length(AllNeurons)
    if iN == 16
        debug = 1;
    end
  NeuronNumber = AllNeurons(iN);
  [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
  disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
  TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, 'DT', [], doitsquare);
    
  filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
  Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
  Expt = Neuron.Expt;
  fileNames{iN} = filename; 
  
    
  dXvalues = unique([Expt.Trials(([Expt.Trials.ce]~=0) & ([Expt.Trials.st]~=0)).dx]);
  dYvalues = unique([Expt.Trials(([Expt.Trials.ce]~=0) & ([Expt.Trials.st]~=0)).dy]);
  if (sum(iN==[16 17 18])>0)
      dXvalues = dXvalues(dXvalues~=0);
  end
  
  SpikeCounts = zeros(length([Expt.Trials]),1);

  for tr = 1: length([Expt.Trials]), 
    if doitsquare
        SpikeCounts(tr) = sqrt(sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime));
    else
        SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
    end
  end
  
  uncorrSpikes = []; blankSpikes = [];
  heatmapCell = cell(length(dXvalues), length(dYvalues)); %zeros(length(dXvalues), length(dYvalues));
  for tr = 1: length([Expt.Trials]), 
      if( (Expt.Trials(tr).ce==0)) % uncorr
          uncorrSpikes = [uncorrSpikes, SpikeCounts(tr)];
      else if( (Expt.Trials(tr).st==0)) % blank
             blankSpikes = [blankSpikes, SpikeCounts(tr)];
          else 
              if(~isempty(find(dXvalues==Expt.Trials(tr).dx)) && ~isempty(find(dYvalues==Expt.Trials(tr).dy)))
                heatmapCell{find(dXvalues==Expt.Trials(tr).dx), find(dYvalues==Expt.Trials(tr).dy)} = [heatmapCell{find(dXvalues==Expt.Trials(tr).dx), find(dYvalues==Expt.Trials(tr).dy)}, SpikeCounts(tr)]  ;
              end
          end
      end
  end
  heatmap = zeros(length(dXvalues), length(dYvalues));
  for i = 1:length(dXvalues)
      for j = 1:length(dYvalues)
        heatmap(i,j) = mean(heatmapCell{i,j});
      end
  end
  disp([num2str(uncorrSpikes) ,' =>' num2str(mean(uncorrSpikes))]), disp([num2str(blankSpikes) ' => ' num2str(mean(blankSpikes))]);
  
  h = 678 + iN;
  figure(h), clf, hold on, set(gcf,'PaperPositionMode','auto'), 
  subplot(2,2,1)
  pcolor(heatmap' - mean(uncorrSpikes))
  shading flat %interp %faceted %flat 
  colorbar
  set(gca, 'XTickLabel', round(dXvalues*10)/10);
  set(gca, 'YTickLabel', round(dYvalues*10)/10);
  
  
  %%
  paddedheatmap = ones(gfs, gfs) .* mean(uncorrSpikes);
  paddedheatmap(1:size(heatmap,1),1:size(heatmap,2)) = heatmap;
  subplot(2,2,4), 
  %pcolor(conv2(Normalheatmap,gfilter))
  MUheatmap = heatmap-mean(uncorrSpikes);
  MUheatmap = (MUheatmap - mean(MUheatmap(:))) ./ range(MUheatmap(:));
  fspace = abs(fftshift(fft2(MUheatmap,gfs,gfs))); %paddedheatmap,gfs,gfs))));
  fspacepeak = max(fspace(:));
  [ax, by] = find(fspace==fspacepeak);
  pcolor(fspace');
  shading flat %interp %faceted %flat 
  colorbar
  hold on
  hL1 = line([size(fspace,1)/2, ax(1)], [size(fspace,2)/2, by(1)], 'Color', [0.0 0.0 0.0]);
  set(hL1,'LineWidth',8);
  hL2 = line([size(fspace,1)/2, ax(1)], [size(fspace,2)/2, by(1)], 'Color', [1.0 1.0 0.0]);
  set(hL1,'LineWidth',4);
  xlim([size(fspace,1)/5 size(fspace,1)*4/5]);
  ylim([size(fspace,2)/5 size(fspace,2)*4/5]);
  
  
  %%
  for x = 1:size(fspace,1)
    for y = 1:size(fspace,2)
      xC = x - size(fspace,1) /2;
      yC = y - size(fspace,2) /2;
      [theta, rho] = cart2pol(xC, yC);
      cartfspace(round((theta * 180 / pi) + 180)+1, round(rho)+1) = fspace(x,y);
    end
  end
  subplot(2,2,2), hold on
  pcolor(cartfspace');
  shading flat
  colorbar
  xlim([1 size(cartfspace,1)]);
  ylim([1 size(cartfspace,2)]);
  
  %%
  URheatmap = heatmap  - mean(uncorrSpikes);
  Normalheatmap = (URheatmap - mean(URheatmap(:))) ./ range(URheatmap(:));
  subplot(2,2,3)
  pcolor(conv2(Normalheatmap',gfilter, 'same'))
  shading flat %interp %faceted %flat 
  colorbar
  speak = max(Normalheatmap(:));
  strough = min(Normalheatmap(:));
  [apx, bpy] = find(Normalheatmap==speak);
  [atx, bty] = find(Normalheatmap==strough);
  % translate then 90 degress rotation then translate back
  apxR = -(bpy-size(Normalheatmap,2)/2) + size(Normalheatmap,2)/2;
  bpyR =  (apx-size(Normalheatmap,1)/2) + size(Normalheatmap,1)/2;
  atxR = -(bty-size(Normalheatmap,2)/2) + size(Normalheatmap,2)/2;
  btyR =  (atx-size(Normalheatmap,1)/2) + size(Normalheatmap,1)/2;
  
  hold on
  hL1 = line([atx(1), apx(1)], [bty(1), bpy(1)], 'Color', [0.0 0.0 0.0]);
  set(hL1,'LineWidth',4);
  hL2 = line([atx(1), apx(1)], [bty(1), bpy(1)], 'Color', [1.0 1.0 0.0]);
  set(hL1,'LineWidth',2);
  hL1 = line([atxR(1), apxR(1)], [btyR(1), bpyR(1)], 'Color', [0.0 0.0 0.0]);
  set(hL1,'LineWidth',8);
  hL2 = line([atxR(1), apxR(1)], [btyR(1), bpyR(1)], 'Color', [0.0 1.0 0.0]);
  set(hL1,'LineWidth',4);
  
  %%
  DIRfilename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.OT.mat');
  if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , DIRfilename),'file'))
    DIRNeuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,DIRfilename));
    DIRExpt = DIRNeuron.Expt;
    DIRvalues = unique([DIRExpt.Trials.('or')]);
  

    DIRSpikeCounts = zeros(length([DIRExpt.Trials]),1);

    for tr = 1: length([DIRExpt.Trials]), 
      if doitsquare
        DIRSpikeCounts(tr) = sqrt(sum([DIRExpt.Trials(tr).Spikes]>=StartTime & [DIRExpt.Trials(tr).Spikes]<=FinishTime));
      else
        DIRSpikeCounts(tr) = sum([DIRExpt.Trials(tr).Spikes]>=StartTime & [DIRExpt.Trials(tr).Spikes]<=FinishTime);
      end
    end
  
    DIRmapCell = cell(length(DIRvalues),1);
    for tr = 1: length([DIRExpt.Trials]), 
      DIRmapCell{find(DIRvalues==DIRExpt.Trials(tr).or)} = [DIRmapCell{find(DIRvalues==DIRExpt.Trials(tr).or)} , DIRSpikeCounts(tr)];
    end
    DIRmap = zeros(length(DIRvalues),1);
    for i = 1: length(DIRvalues), 
      DIRmap(i) = mean(DIRmapCell{i});
    end
     
%    subplot(2,2,2)
%    hp = polar([DIRvalues(end) DIRvalues] .* pi ./ 180, [DIRmap(end); DIRmap]');
%    set(hp, 'LineWidth',4),
  end
  
  %%
  ORfilename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rls.OT.mat');
  if(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , ORfilename),'file'))
    ORNeuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,ORfilename));
    ORExpt = ORNeuron.Expt;
    if (isfield(ORExpt.Trials,'or')) 
      ORvalues = unique([ORExpt.Trials.('or')]);
  

      ORSpikeCounts = zeros(length([ORExpt.Trials]),1);

      for tr = 1: length([ORExpt.Trials]), 
        if doitsquare
          ORSpikeCounts(tr) = sqrt(sum([ORExpt.Trials(tr).Spikes]>=StartTime & [ORExpt.Trials(tr).Spikes]<=FinishTime));
        else
          ORSpikeCounts(tr) = sum([ORExpt.Trials(tr).Spikes]>=StartTime & [ORExpt.Trials(tr).Spikes]<=FinishTime);
        end
      end
  
      ORmapCell = cell(length(ORvalues),1);
      for tr = 1: length([ORExpt.Trials]), 
        ORmapCell{find(ORvalues==ORExpt.Trials(tr).or)} = [ORmapCell{find(ORvalues==ORExpt.Trials(tr).or)} , ORSpikeCounts(tr)];
      end
      ORmap = zeros(length(ORvalues),1);
      for i = 1: length(ORvalues), 
        ORmap(i) = mean(ORmapCell{i});
      end
     
 %     subplot(2,2,2), hold on
 %     hp = polar([ORvalues ORvalues+180] .* pi ./ 180, [ORmap; ORmap]', 'r');
 %     set(hp, 'LineWidth',4),
    end
    
  end
  
  
  %%
  ebX=[]; cbX = [];
  ebY=[]; cbY = [];
  for i = 1:length(dXvalues)
        ebX(i) = mean(SpikeCounts([Expt.Trials.dx]==dXvalues(i)));
        cbX(i) = sum([Expt.Trials.dx]==dXvalues(i));
  end
  for i = 1:length(dYvalues)
        ebY(i) = mean(SpikeCounts([Expt.Trials.dy]==dYvalues(i)));
        cbY(i) = sum([Expt.Trials.dy]==dYvalues(i));
  end

  tisX = [];
  for i = 1:floor(length(ebX)/2)
    a = i;
    b = 1 + length(ebX) - i;
    tisX(i) = (sum(ebX(a) .* cbX(a)) ./ mean(cbX(a)) - sum(ebX(b) .* cbX(b)) ./ mean(cbX(b)) ) / ...
     (sum(ebX(a) .* cbX(a)) ./ mean(cbX(a)) + sum(ebX(b) .* cbX(b)) ./ mean(cbX(b)) );

    b = 1 + length(ebY) - i;
    tisY(i) = (sum(ebY(a) .* cbY(a)) ./ mean(cbY(a)) - sum(ebY(b) .* cbY(b)) ./ mean(cbY(b)) ) / ...
     (sum(ebY(a) .* cbY(a)) ./ mean(cbY(a)) + sum(ebY(b) .* cbY(b)) ./ mean(cbY(b)) );
  end

TIX = max(tisX);
if (TIX < abs(min(tisX)))
    TIX = min(tisX);
end

TIY = max(tisY);
if (TIY < abs(min(tisY)))
    TIY = min(tisY);
end

if(isempty(TIX))
    TIX = -998;
end
if(isempty(TIY))
    TIY = -998;
end

TIXY(iN,:) = [TIX, TIY]; 
set(h, 'Units', 'pixels');
pW = 440; pH = 330;
a = (mod(iN,6)) * pW;
b = (fix(iN/6)) * pH;
set(h, 'OuterPosition', [0+a 1280-b pW pH])
print(h, '-dpsc', '-r150', '-zbuffer', ['../figs/DXDY', num2str(h), '.eps']);
end

%%
figure(18), clf, hold on,
scatter(TIXY(:,1), TIXY(:,2), 'filled')
scatter(TIXY(:,1), TI, 'filled')
refline(0)
refline(1)
reflinexy(0,1);
