
clear, clc
DataPath = '/bgc/data/';

% % % % DID
load('../AllDIDNeurons.mat');
AllNeurons = AllDIDNeurons;
clear AllDIDNeurons;
FileType = 'DID';
StimulusType = 'cylinder';


% % % TWO
% load('../AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWOPsychDays;
% FileType = 'TWO';
% StimulusType = 'cylinder';

%par
for iN= 1:length(AllNeurons), 
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , num2str(rem(now,1))));
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'); 
    filenames{iN} = strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    
    if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', filename), 'file')==2)
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;

        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end

        switch FileType
            case 'TWO'
                cip = ([Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive) ;
                cip0= ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
                cin = ([Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                cin0= ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                wip = ([Expt.Trials(:).dx]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                wip0= ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                win = ([Expt.Trials(:).dx]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
                win0= ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
            case 'DID'
                IdPref = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0);
                IdNull = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0);
                
                RespPs = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
                RespNg = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                
                RespDP = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]== 1);
                RespDN = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==-1);
        end
        lv = zeros(length(Expt.Trials), 2000); lh = zeros(length(Expt.Trials), 2000);
        %si = round(sum(Expt.Header.emtimes<0)/10)+1; % stimulus onset
        %disp([' ----------------------------- ', num2str(si)]);
        for i = 1: length(Expt.Trials)
            validpoints = Expt.Header.emtimes>0 & Expt.Header.emtimes < 25000;
            if strcmpi(MonkeyName, 'icarus')
%                lh(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),1);
%                lv(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),3);
                 lh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,1);
                 lv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,3);
            else
                %lh(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),2);
                %lv(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),4);
                lh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,2);
                lv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,4);
            end
            if (~isempty([Expt.Trials(i).Saccades]))
                saccads(i,:,:) = [mean([Expt.Trials(i).Saccades.dir]); mean([Expt.Trials(i).Saccades.size])]; 
            else
                saccads(i,:,:) = [0; 0];
            end
        end
        figure(1892), hist(saccads([Expt.Trials(:).RespDir]~=0,1));
        switch(FileType)
                case 'TWO'
               
                case 'DID'
                    IdPv = lv(IdPref, 1:2000);
                    IdPh = lh(IdPref, 1:2000);
                    IdNv = lv(IdNull, 1:2000);
                    IdNh = lh(IdNull, 1:2000);

                    RPsv = lv(RespPs, 1:2000);
                    RPsh = lh(RespPs, 1:2000);
                    RNsv = lv(RespNg, 1:2000);
                    RNsh = lh(RespNg, 1:2000); 
                    
                    RePv = lv(RespDP, 1:2000);
                    RePh = lh(RespDP, 1:2000);
                    ReNv = lv(RespDN, 1:2000);
                    ReNh = lh(RespDN, 1:2000); 
                    
                    % BASELINE Correction
                    IdPv = IdPv - mean(mean(IdPv(:,1:1000)));
                    IdPh = IdPh - mean(mean(IdPh(:,1:1000)));
                    IdNv = IdNv - mean(mean(IdNv(:,1:1000)));
                    IdNh = IdNh - mean(mean(IdNh(:,1:1000)));
                    
                    RPsv = RPsv - mean(mean(RPsv(:,1:1000)));
                    RPsh = RPsh - mean(mean(RPsh(:,1:1000)));
                    RNsv = RNsv - mean(mean(RNsv(:,1:1000)));
                    RNsh = RNsh - mean(mean(RNsh(:,1:1000)));
                    
                    RePv = RePv - mean(mean(RePv(:,1:1000)));
                    RePh = RePh - mean(mean(RePh(:,1:1000)));
                    ReNv = ReNv - mean(mean(ReNv(:,1:1000)));
                    ReNh = ReNh - mean(mean(ReNh(:,1:1000)));
                    
                    theta = Expt.Stimvals.or*pi / 180.0;
                    rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                    
                    r = [mean(IdPh); mean(IdPv)]' * rotationMatrix;
                    rIdPh = r(:,1)'; rIdPv = r(:,2)';
                    
                    r = [mean(IdNh); mean(IdNv)]' * rotationMatrix;
                    rIdNh = r(:,1)'; rIdNv = r(:,2)';

                    r = [mean(RPsh); mean(RPsv)]' * rotationMatrix;
                    rRPsh = r(:,1)'; rRPsv = r(:,2)';

                    r = [mean(RNsh); mean(RNsv)]' * rotationMatrix;
                    rRNsh = r(:,1)'; rRNsv = r(:,2)';

                    r = [mean(RePh); mean(RePv)]' * rotationMatrix;
                    rRePh = r(:,1)'; rRePv = r(:,2)';

                    r = [mean(ReNh); mean(ReNv)]' * rotationMatrix;
                    rReNh = r(:,1)'; rReNv = r(:,2)';

                    % a = rIdPv - rIdNv; b = rIdPh - rIdNh;
                    a = rIdPh; b = rIdNh;
                    if a(1600) > b(1600)
                        crossId = max(find(a(1:1600) < b(1:1600)));
                    else
                         crossId = max(find(a(1:1600) > b(1:1600)));
                    end
                    
                    a = rRPsh; b = rRNsh;
                    if a(1600) > b(1600)
                        crossCh = max(find(a(1:1600) < b(1:1600)));
                    else
                         crossCh = max(find(a(1:1600) > b(1:1600)));
                    end
                    
                    a = rRePh; b = rReNh;
                    if a(1600) > b(1600)
                        crossCe = max(find(a(1:1600) < b(1:1600)));
                    else
                         crossCe = max(find(a(1:1600) > b(1:1600)));
                    end
                    
                    scaling(iN) = Expt.Header.emtimes(end) / length(Expt.Header.emtimes) / 10; %Expt.Header.emtimes(end);
                    disp([ ' -----------------------------' , num2str(scaling(iN))]);
                    % =================   Micros Saccades
                    SacP  = saccads(IdPref,:);
                    SacN  = saccads(IdNull,:);
                    SacRP = saccads(RespPs,:);
                    SacRN = saccads(RespNg,:);
                    SacReP= saccads(RespDP,:);
                    SacReN= saccads(RespDN,:);
                    
                    EyeTrack{iN} = {IdPv, IdPh, IdNv, IdNh, RPsv, RPsh, RNsv, RNsh, ...
                        rIdPv, rIdPh, rIdNv, rIdNh, rRPsv, rRPsh, rRNsv, rRNsh, ...
                        crossId, crossCh, scaling(iN), ...
                        Expt.Stimvals.or, MonkeyName, Expt.Stimvals.xo, Expt.Stimvals.yo, ... 
                        median([Expt.Trials(:).sz]), Expt.Stimvals.bo, SacP, SacN, SacRP, SacRN, ...
                        rRPsv, rRPsh, rRNsv, rRNsh, crossCe, SacReP, SacReN};
        end
    else
        disp('Wow!')
    end
end

%% Graphics

for iN= 1:length(AllNeurons),
   if (0<size(EyeTrack{iN},1))
        SaccPC(iN) = size(EyeTrack{iN}{26},1);
        SaccNC(iN) = size(EyeTrack{iN}{27},1);
        SacRPC(iN) = size(EyeTrack{iN}{28},1);
        SacRNC(iN) = size(EyeTrack{iN}{29},1);
        SacRePC(iN)= size(EyeTrack{iN}{35},1);
        SacReNC(iN)= size(EyeTrack{iN}{36},1);

        SaccP(iN,:) = mean(EyeTrack{iN}{26});
        SaccN(iN,:) = mean(EyeTrack{iN}{27});
        SacRP(iN,:) = mean(EyeTrack{iN}{28});
        SacRN(iN,:) = mean(EyeTrack{iN}{29});
        SacReP(iN,:)= mean(EyeTrack{iN}{35});
        SacReN(iN,:)= mean(EyeTrack{iN}{36});

   end
end
%%
figure(368), clf, 
scatter(SaccPC, SaccNC, 'filled')
refline(1);

figure(378), clf, 
scatter(SacRPC, SacRNC, 'filled')
refline(1);

figure(388), clf, 
scatter(SacRePC, SacReNC, 'filled')
refline(1);

figure(318), clf, 
scatter(SaccP(:,1), SaccN(:,1), 'filled')
refline(1);

figure(328), clf, 
scatter(SacRP(:,1), SacRN(:,1), 'filled')
refline(1);

figure(338), clf, 
scatter(SaccP(:,2), SaccN(:,2), 'filled')
refline(1);

figure(348), clf, 
scatter(SacRP(:,2), SacRN(:,2), 'filled')
refline(1);



%%
figure(1112), clf, hold on,
m = []; n = [];
for iN= 1:length(AllNeurons),
    if ~isempty(EyeTrack{iN})
        scale100(iN) = EyeTrack{iN}{19};
        
        m(iN) = EyeTrack{iN}{17} * scale100(iN);
        n(iN) = EyeTrack{iN}{18} * scale100(iN);
        
        mt(iN, :) = EyeTrack{iN}{10};
        nt(iN, :) = EyeTrack{iN}{12};
        
        mtc(iN, :) = EyeTrack{iN}{14};
        ntc(iN, :) = EyeTrack{iN}{16};
    else
        m(iN) = -1;
        n(iN) = -1;
        %mt(iN,:)= 0;
        %nt(iN,:)= 0;
    end
end

Xaxis = [1:2000] * mean(scale100(scale100>0));

clickscatter(m, n, 2, [], filenames);
refline(1);

figure(1212), clf, hold on,
plot(Xaxis, mean(mt), 'r', 'LineWidth', 2), 
plot(Xaxis, mean(nt), 'b', 'LineWidth', 2);
plot(Xaxis, mean(mtc), 'm', 'LineWidth', 2), 
plot(Xaxis, mean(ntc), 'k', 'LineWidth', 2);

set(gca, 'XGrid', 'on');
%xlim([0 2000]);
xtl = [0, 500, 1000, 1500, 2000, 2500];
set(gca, 'XTick', xtl);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);





%%
for iN= 1:length(AllNeurons),
        figure(11112), clf, hold on,
%         plot(mean(EyeTrack{iN}{1}) - mean(EyeTrack{iN}{3}), 'r.', 'LineWidth', 2);
%         plot(mean(EyeTrack{iN}{2}) - mean(EyeTrack{iN}{4}), 'g.', 'LineWidth', 2);
%         plot(mean(EyeTrack{iN}{5}) - mean(EyeTrack{iN}{7}), 'b.', 'LineWidth', 2);
%         plot(mean(EyeTrack{iN}{6}) - mean(EyeTrack{iN}{8}), 'm.', 'LineWidth', 2);

        %plot(EyeTrack{iN}{9}  - EyeTrack{iN}{11}, 'r', 'LineWidth', 2);
        %plot(EyeTrack{iN}{10} - EyeTrack{iN}{12}, 'g', 'LineWidth', 2);
%        plot(EyeTrack{iN}{13} - EyeTrack{iN}{14}, 'b', 'LineWidth', 2);
%        plot(EyeTrack{iN}{15} - EyeTrack{iN}{16}, 'm', 'LineWidth', 2);
        
        %plot(EyeTrack{iN}{9},  'r.', 'LineWidth', 1);
        plot(Xaxis, EyeTrack{iN}{10}, 'g.', 'LineWidth', 1);
        %plot(EyeTrack{iN}{11}, 'b--', 'LineWidth', 1);
        plot(Xaxis, EyeTrack{iN}{12}, 'k.', 'LineWidth', 1);
        reflinexy(EyeTrack{iN}{17} * mean(scale100(scale100>0)), 1);
        
        plot(Xaxis, EyeTrack{iN}{14}, 'r--', 'LineWidth', 1);
        plot(Xaxis, EyeTrack{iN}{16}, 'b--', 'LineWidth', 1);
        reflinexy(EyeTrack{iN}{18} * mean(scale100(scale100>0)), 1, 'LineWidth', 3);
        
        pause
end


%% 
for iN= 1:length(AllNeurons),
        figure(11113), clf, hold on,
        plot(Xaxis, EyeTrack{iN}{9},  'r', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{10}, 'g', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{11}, 'r--', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{12}, 'g--', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{13}, 'c', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{14}, 'm', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{15}, 'c--', 'LineWidth', 2);
        plot(Xaxis, EyeTrack{iN}{16}, 'm--', 'LineWidth', 2);
        pause
end

