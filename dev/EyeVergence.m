clear, clc

DataPath = GetDataPath('server');

% % % % DID
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('DID');

% experiments with valid eyetracking for  both eyes
BinocularExperiments = {'dae014';'dae015';'dae018';'dae054';'dae149';'dae151';'dae447';'dae4472';'dae449';'dae4492';};

%load ~/Desktop/matlab.mat EyeTrack

%par
for iN= 1:length(AllNeurons), 
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , num2str(rem(now,1))));
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'); 
    %disp(filename);
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
        
        pD(iN) = PreferredCylinderRotationDirection(MonkeyName, NeuronNumber, ClusterName, FileType, 0);
        
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
                if (pD == 2)
                    IdPref = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0);
                    IdNull = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0);
                    RespPs = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                    RespNg = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
                else
                    IdPref = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir]~=0);
                    IdNull = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir]~=0);
                    RespPs = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToPositive);
                    RespNg = ([Expt.Trials(:).dx]==0 & [Expt.Trials(:).RespDir]==ResponseToNegative);
                end
        end
        lv = zeros(length(Expt.Trials), 2000); lh = zeros(length(Expt.Trials), 2000);
        rv = zeros(length(Expt.Trials), 2000); rh = zeros(length(Expt.Trials), 2000);
        for i = 1: length(Expt.Trials)
            validpoints = Expt.Header.emtimes>0 & Expt.Header.emtimes < 25000;
            lh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,1);
            lv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,3);
            rh(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,2);
            rv(i, 1:sum(validpoints)) = Expt.Trials(i).EyeData(validpoints,4);
        end
        
% %         % sanity! 
%         figure(4321), clf, hold on
%         plot(mean(lh), 'r', 'LineWidth', 2)
%         plot(mean(rh), 'b', 'LineWidth', 2)
%         plot(mean(lv), 'm', 'LineWidth', 2)
%         plot(mean(rv), 'k', 'LineWidth', 2)
%         
%         figure(5432), clf, hold on, 
%         plot(mean(lh) - mean(rh), 'r--', 'LineWidth', 2);
%         plot(mean(lv) - mean(rv), 'b--', 'LineWidth', 2);
%          
        switch(FileType)
                case 'TWO'
               
                case 'DID'
                    IdPlv = lv(IdPref, 1:2000);
                    IdPlh = lh(IdPref, 1:2000);
                    IdNlv = lv(IdNull, 1:2000);
                    IdNlh = lh(IdNull, 1:2000);

                    RPslv = lv(RespPs, 1:2000);
                    RPslh = lh(RespPs, 1:2000);
                    RNslv = lv(RespNg, 1:2000);
                    RNslh = lh(RespNg, 1:2000); 

                    IdPrv = rv(IdPref, 1:2000);
                    IdPrh = rh(IdPref, 1:2000);
                    IdNrv = rv(IdNull, 1:2000);
                    IdNrh = rh(IdNull, 1:2000);

                    RPsrv = rv(RespPs, 1:2000);
                    RPsrh = rh(RespPs, 1:2000);
                    RNsrv = rv(RespNg, 1:2000);
                    RNsrh = rh(RespNg, 1:2000); 

%                     % BASELINE Correction
                    IdPlv = IdPlv - mean(mean(IdPlv(:,1:100)));
                    IdPlh = IdPlh - mean(mean(IdPlh(:,1:100)));
                    IdNlv = IdNlv - mean(mean(IdNlv(:,1:100)));
                    IdNlh = IdNlh - mean(mean(IdNlh(:,1:100)));
                    
                    RPslv = RPslv - mean(mean(RPslv(:,1:100)));
                    RPslh = RPslh - mean(mean(RPslh(:,1:100)));
                    RNslv = RNslv - mean(mean(RNslv(:,1:100)));
                    RNslh = RNslh - mean(mean(RNslh(:,1:100)));
                    
                    IdPrv = IdPrv - mean(mean(IdPrv(:,1:100)));
                    IdPrh = IdPrh - mean(mean(IdPrh(:,1:100)));
                    IdNrv = IdNrv - mean(mean(IdNrv(:,1:100)));
                    IdNrh = IdNrh - mean(mean(IdNrh(:,1:100)));
                    
                    RPsrv = RPsrv - mean(mean(RPsrv(:,1:100)));
                    RPsrh = RPsrh - mean(mean(RPsrh(:,1:100)));
                    RNsrv = RNsrv - mean(mean(RNsrv(:,1:100)));
                    RNsrh = RNsrh - mean(mean(RNsrh(:,1:100)));

                    scaling(iN) = Expt.Header.emtimes(end) / length(Expt.Header.emtimes) / 10; %Expt.Header.emtimes(end);
                    disp([ ' -----------------------------' , num2str(scaling(iN))]);
                                        
                    EyeTrack{iN} = {
                        IdPlv, IdPlh, IdNlv, IdNlh, RPslv, RPslh, RNslv, RNslh, ...
                        IdPrv, IdPrh, IdNrv, IdNrh, RPsrv, RPsrh, RNsrv, RNsrh, ...
                        scaling(iN), Expt.Stimvals.or, MonkeyName, Expt.Stimvals.xo, Expt.Stimvals.yo, ... 
                        median([Expt.Trials(:).sz]), Expt.Stimvals.bo};
        end
    else
        disp('Wow!')
    end
end

%% Graphics

iNc = 0;
for iN= 1:length(AllNeurons),
    if ~isempty(EyeTrack{iN}) && (sum(strcmpi(AllNeurons{iN}, BinocularExperiments(:)))) 
            iNc = iNc + 1;
            scale100(iNc) = EyeTrack{iN}{17};
        
            vIdP(iNc, :) = mean(EyeTrack{iN}{2} - EyeTrack{iN}{10});
            vIdN(iNc, :) = mean(EyeTrack{iN}{4} - EyeTrack{iN}{12});

            vRP(iNc, :) = mean(EyeTrack{iN}{6} - EyeTrack{iN}{14});
            vRN(iNc, :) = mean(EyeTrack{iN}{8} - EyeTrack{iN}{16});
    
            vvIdP(iNc, :) = mean(EyeTrack{iN}{1} - EyeTrack{iN}{9});
            vvIdN(iNc, :) = mean(EyeTrack{iN}{3} - EyeTrack{iN}{11});

            vvRP(iNc, :) = mean(EyeTrack{iN}{5} - EyeTrack{iN}{13});
            vvRN(iNc, :) = mean(EyeTrack{iN}{7} - EyeTrack{iN}{15});
            
            disp(pD(iN));
            roc(iNc) = rocs(iN, 4);
    end
end

Xaxis = [1:2000] * mean(scale100(scale100>0));

figure(2212), clf, hold on,
plot(Xaxis, mean(vIdP), 'r', 'LineWidth', 2), 
plot(Xaxis, mean(vIdN), 'b', 'LineWidth', 2);
plot(Xaxis, mean(vRP), 'm', 'LineWidth', 2), 
plot(Xaxis, mean(vRN), 'k', 'LineWidth', 2);

set(gca, 'XGrid', 'on');
%xlim([0 2000]);
xtl = [0, 500, 1000, 1500, 2000, 2500];
set(gca, 'XTick', xtl);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);
refline(0)

figure(2112), clf, hold on,
plot(Xaxis, mean(vIdP - vIdN), 'r', 'LineWidth', 2), 
plot(Xaxis, mean(vRP - vRN), 'b', 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([0 2500]);
xtl = [0, 500, 1000, 1500, 2000, 2500];
set(gca, 'XTick', xtl);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);
ylim([-0.20002 0.700007]);
refline(0);

[a, b, c, d] = ttest2(mean(vIdP(:,1:1350)),mean(vIdN(:,1:1350)))
[a, b, c, d] = ttest2(mean(vRP(:,1:1350)),mean(vRN(:,1:1350)))

figure(2113), clf, hold on,
plot(Xaxis, mean(vvIdP - vvIdN), 'r', 'LineWidth', 2), 
plot(Xaxis, mean(vvRP - vvRN), 'b', 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([0 2500]);
xtl = [0, 500, 1000, 1500, 2000, 2500];
set(gca, 'XTick', xtl);
set(gca, 'XTickLabel', {num2str(xtl')});
title(FileType);
refline(0)

[a, b, c, d] = ttest2(mean(vvIdP(:,1:1350)),mean(vvIdN(:,1:1350)))
[a, b, c, d] = ttest2(mean(vvRP(:,1:1350)),mean(vvRN(:,1:1350)))


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

