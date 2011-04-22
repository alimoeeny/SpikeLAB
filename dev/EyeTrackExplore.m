% Eye track explore
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
        end
        lv = zeros(length(Expt.Trials), 2000); lh = zeros(length(Expt.Trials), 2000);
        lvSacc = zeros(length(Expt.Trials), 2000); lhSacc = zeros(length(Expt.Trials), 2000);

        %disp(Expt.Stimvals.wr);
        saccads = {}; signal = []; SacPool = {}; AllSacc = {};
        si = round(sum(Expt.Header.emtimes<0)/10)+1; % /10 added 04/04/2011
        for i = 1: length(Expt.Trials)
            signal(i) = Expt.Trials(i).dx;
            if ~isempty(Expt.Trials(i).EyeData) && Expt.Trials(i).RespDir ~=0
                done = 0;
                for tid = 1: length(Expt.Trials(i).Saccades),
                    AllSacc{i}(tid) = Expt.Trials(i).Saccades(tid);
                    if(done == 0 && Expt.Trials(i).Saccades(tid).start<25000 && Expt.Trials(i).Saccades(tid).start>20000 && Expt.Trials(i).Saccades(tid).size >0.5)
                        SacPool{end+1} = Expt.Trials(i).Saccades(tid);
                        SacPool{end}.TrialIdx = i; 
                        done = 1;

                        if strcmpi(MonkeyName, 'icarus')
                            lh(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),1);
                            lv(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),3);

                            saccads{i} = Expt.Trials(i).Saccades(tid); 
                            saccI = round(saccads{i}.start / 10);
                            saccI = find(Expt.Header.emtimes/10>saccI,1);
                            lhSacc(
                            
                            
                            
                            This file does not work since it assumes a wrong sample rate for eye data
                            emtimes contains the actual time stamps and is used incorrectly here
                            refere to EyeMovement.m for and updated version
                            
                            
                            i,1:1+length(Expt.Trials(i).EyeData)-saccI) = Expt.Trials(i).EyeData((saccI:end),1);
                            lvSacc(i,1:1+length(Expt.Trials(i).EyeData)-saccI) = Expt.Trials(i).EyeData((saccI:end),3);
                        else
                            lh(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),2);
                            lv(i,1:1+size(Expt.Trials(i).EyeData,1)-si) = Expt.Trials(i).EyeData((si:end),4);

                            saccads{i} = Expt.Trials(i).Saccades(tid); 
                            saccI = round(saccads{i}.start / 10);
                            saccI = find(Expt.Header.emtimes/10>saccI,1);
                            lhSacc(i,1:1+length(Expt.Trials(i).EyeData)-saccI) = Expt.Trials(i).EyeData((saccI:end),2);
                            lvSacc(i,1:1+length(Expt.Trials(i).EyeData)-saccI) = Expt.Trials(i).EyeData((saccI:end),4);
                        end
                        SacPool{end}.lh = lh(i,:);
                        SacPool{end}.lv = lv(i,:);
                    end
                end
            else
                lv(i,:) = 0;
                lh(i,:) = 0;
            end
        end

        switch(FileType)
            case 'TWO'
                cpt = lv(cip,1:1800);
                cpth= lh(cip,1:1800);
                wnt = lv(win,1:1800);
                wnth= lh(win,1:1800);
                cnt = lv(cin,1:1800);
                cnth= lh(cin,1:1800);
                wpt = lv(wip,1:1800);
                wpth= lh(wip,1:1800);

                cpt0 = lv(cip0,1:1800);
                cpth0= lh(cip0,1:1800);
                wnt0 = lv(win0,1:1800);
                wnth0= lh(win0,1:1800);
                cnt0 = lv(cin0,1:1800);
                cnth0= lh(cin0,1:1800);
                wpt0 = lv(wip0,1:1800);
                wpth0= lh(wip0,1:1800);

                cptSacc = lvSacc(cip,1:2000);
                cpthSacc= lhSacc(cip,1:2000);
                wntSacc = lvSacc(win,1:2000);
                wnthSacc= lhSacc(win,1:2000);
                cntSacc = lvSacc(cin,1:2000);
                cnthSacc= lhSacc(cin,1:2000);
                wptSacc = lvSacc(wip,1:2000);
                wpthSacc= lhSacc(wip,1:2000);

                cptSacc0 = lvSacc(cip0,1:2000);
                cpthSacc0= lhSacc(cip0,1:2000);
                wntSacc0 = lvSacc(win0,1:2000);
                wnthSacc0= lhSacc(win0,1:2000);
                cntSacc0 = lvSacc(cin0,1:2000);
                cnthSacc0= lhSacc(cin0,1:2000);
                wptSacc0 = lvSacc(wip0,1:2000);
                wpthSacc0= lhSacc(wip0,1:2000);

                EyeTrack{iN} = {cpt, cpth, wnt, wnth, cnt, cnth, wpt, wpth, Expt.Stimvals.or, MonkeyName, Expt.Stimvals.xo, ... %11
                            Expt.Stimvals.yo, median([Expt.Trials(:).sz]), Expt.Stimvals.bo, cptSacc, cpthSacc, wntSacc,... %17
                            wnthSacc, cntSacc, cnthSacc, wptSacc, wpthSacc, saccads, signal, SacPool, cip, cin, wip, win, ... %29
                            cip0, cin0, wip0, win0, cpt0, cpth0, wnt0, wnth0, cnt0, cnth0, wpt0, wpth0, cptSacc0, cpthSacc0, ... %43
                            wntSacc0, wnthSacc0, cntSacc0, cnthSacc0, wptSacc0, wpthSacc0, AllSacc};
            case 'DID'
                IdPv = lv(IdPref, 1:1800);
                IdPh = lh(IdPref, 1:1800);
                IdNv = lv(IdNull, 1:1800);
                IdNh = lh(IdNull, 1:1800);

                RPsv = lv(RespPs, 1:1800);
                RPsh = lh(RespPs, 1:1800);
                RNsv = lv(RespNg, 1:1800);
                RNsh = lh(RespNg, 1:1800);

                EyeTrack{iN} = {IdPv, IdPh, IdNv, IdNh, RPsv, RPsh, RNsv, RNsh, Expt.Stimvals.or, MonkeyName, Expt.Stimvals.xo, ...
                    Expt.Stimvals.yo, median([Expt.Trials(:).sz]), Expt.Stimvals.bo, AllSacc};

        end
    end
    
    %figure(832), clf, hold on,
    %for i = 1:length(ecipd)
    %    plot(eExpt.Trials(ecipd(i)).Eyevals.lh)
    %     plot(eExpt.Trials(ecipd(i)).Eyevals.lv, 'r')
    %     cpt(i,1:length(eExpt.Trials(ecipd(i)).Eyevals.lv)) = eExpt.Trials(ecipd(i)).Eyevals.lv;
    %     cpth(i,1:length(eExpt.Trials(ecipd(i)).Eyevals.lh)) = eExpt.Trials(ecipd(i)).Eyevals.lh;
    %    plot(lv(ecipd(i),:), 'r')
    %end
    % 
    % for i = 1:length(ewind)
    % %    plot(eExpt.Trials(ewind(i)).Eyevals.lh)
    %     plot(eExpt.Trials(ewind(i)).Eyevals.lv, 'k')
    %     wnt(i,1:length(eExpt.Trials(ewind(i)).Eyevals.lv)) = eExpt.Trials(ewind(i)).Eyevals.lv;
    %     wnth(i,1:length(eExpt.Trials(ewind(i)).Eyevals.lh)) = eExpt.Trials(ewind(i)).Eyevals.lh;
    % end

    % for i = 1:length(ecind)
    % %    plot(eExpt.Trials(ecipd(i)).Eyevals.lh)
    %     plot(eExpt.Trials(ecind(i)).Eyevals.lv, 'r')
    %     cnt(i,1:length(eExpt.Trials(ecind(i)).Eyevals.lv)) = eExpt.Trials(ecind(i)).Eyevals.lv;
    %     cnth(i,1:length(eExpt.Trials(ecind(i)).Eyevals.lh)) = eExpt.Trials(ecind(i)).Eyevals.lh;
    % end
    % for i = 1:length(ewipd)
    % %    plot(eExpt.Trials(ewind(i)).Eyevals.lh)
    %     plot(eExpt.Trials(ewipd(i)).Eyevals.lv, 'k')
    %     wpt(i,1:length(eExpt.Trials(ewipd(i)).Eyevals.lv)) = eExpt.Trials(ewipd(i)).Eyevals.lv;
    %     wpth(i,1:length(eExpt.Trials(ewipd(i)).Eyevals.lh)) = eExpt.Trials(ewipd(i)).Eyevals.lh;
    % end
    %xlim([0 1800])

%     figure(938), clf, hold on,
%     plot(mean(cpt));
%     plot(mean(wnt),'k');
%     plot(mean(cnt),'b--');
%     plot(mean(wpt),'k--');
% 
%     %figure(939), clf, hold on,
%     plot(mean(cpth),'m', 'LineWidth',2);
%     plot(mean(wnth),'k', 'LineWidth',2);
%     plot(mean(cnth),'m--', 'LineWidth',2);
%     plot(mean(wpth),'k--', 'LineWidth',2);
%     xlim([0 1800])
% 
%     figure(940), clf, hold on,
%     plot(mean(cpt(:,1:1800)) - mean(wnt(:,1:1800)),'r', 'LineWidth',2);
%     plot(mean(cnt(:,1:1800)) - mean(wpt(:,1:1800)),'g', 'LineWidth',2);
%     plot(mean(cpth(:,1:1800)) - mean(wnth(:,1:1800)),'b--', 'LineWidth',2);
%     plot(mean(cnth(:,1:1800)) - mean(wpth(:,1:1800)),'m--', 'LineWidth',2);
%     xlim([0 1500])
%     x = round(2 / eExpt.Header.CRrates(1));
%     % line([x -5],[x 5])


    % plot(eExpt.Trials([eExpt.Trials(:).id]==9162).Eyevals.lh)
    % hold on
    % plot(eExpt.Trials([eExpt.Trials(:).id]==9162).Eyevals.lv, 'r')

end

%% visualization
plotFigures = 0;

for iN = length(EyeTrack):-1:1
     switch(FileType)
         case 'DID'
            IdPv(iN, :) = mean(squeeze(EyeTrack{iN}{1}));
            IdPh(iN, :)  = mean(squeeze(EyeTrack{iN}{2}));
            IdNv(iN, :)  = mean(squeeze(EyeTrack{iN}{3}));
            IdNh(iN, :)  = mean(squeeze(EyeTrack{iN}{4}));
            RPsv(iN, :)  = mean(squeeze(EyeTrack{iN}{5}));
            RPsh(iN, :)  = mean(squeeze(EyeTrack{iN}{6}));
            RNsv(iN, :)  = mean(squeeze(EyeTrack{iN}{7}));
            RNsh(iN, :)  = mean(squeeze(EyeTrack{iN}{8}));

            figure(1100), clf, hold on, 
            plot(squeeze(IdPv(iN, :)), 'r', 'LineWidth', 2);
            plot(squeeze(IdPh(iN, :)), 'b', 'LineWidth', 2);
            plot(squeeze(IdNv(iN, :)), 'r-', 'LineWidth', 2);
            plot(squeeze(IdNh(iN, :)), 'b-', 'LineWidth', 2);
            plot(squeeze(RPsv(iN, :)), 'c.-', 'LineWidth', 2);
            plot(squeeze(RPsh(iN, :)), 'm.-', 'LineWidth', 2);
            plot(squeeze(RNsv(iN, :)), 'b--', 'LineWidth', 2);
            plot(squeeze(RNsh(iN, :)), 'g--', 'LineWidth', 2);
            
            or  = EyeTrack{iN}{9};
            MonkeyName = EyeTrack{iN}{10};
            xo = EyeTrack{iN}{11};
            yo = EyeTrack{iN}{12};
            sz = EyeTrack{iN}{13};
            bo = EyeTrack{iN}{14};
            AllSacc = EyeTrack{iN}{15};
            
         case 'TWO'
            cpt  = EyeTrack{iN}{1};
            cpth = EyeTrack{iN}{2};
            wnt  = EyeTrack{iN}{3};
            wnth = EyeTrack{iN}{4};
            cnt  = EyeTrack{iN}{5};
            cnth = EyeTrack{iN}{6};
            wpt  = EyeTrack{iN}{7};
            wpth = EyeTrack{iN}{8};
            Aor(iN) = EyeTrack{iN}{9};
            Axo(iN) = EyeTrack{iN}{11};
            Ayo(iN) = EyeTrack{iN}{12};
            Asz(iN) = EyeTrack{iN}{13};
            Abo(iN) = EyeTrack{iN}{14};

             cptSacc  = EyeTrack{iN}{15};
             cpthSacc = EyeTrack{iN}{16};
             wntSacc  = EyeTrack{iN}{17};
             wnthSacc = EyeTrack{iN}{18};
             cntSacc  = EyeTrack{iN}{19};
             cnthSacc = EyeTrack{iN}{20};
             wptSacc  = EyeTrack{iN}{21};
             wpthSacc = EyeTrack{iN}{22};

             saccads  = EyeTrack{iN}{23};
             signal   = EyeTrack{iN}{24};
             SacPool  = EyeTrack{iN}{25};

             cip      = EyeTrack{iN}{26};
             cin      = EyeTrack{iN}{27};
             wip      = EyeTrack{iN}{28};
             win      = EyeTrack{iN}{29};

             cip0     = EyeTrack{iN}{30};
             cin0     = EyeTrack{iN}{31};
             wip0     = EyeTrack{iN}{32};
             win0     = EyeTrack{iN}{33};

             AllSacc  = EyeTrack{iN}{50};
             
             %%%% DIFF %%%%%

             cpt = diff(cpt')';
             wpt = diff(wpt')';
             cnt = diff(cnt')';
             wnt = diff(wnt')';
             cpth = diff(cpth')';
             wpth = diff(wpth')';
             cnth = diff(cnth')';
             wnth = diff(wnth')';

             %%%%%%%%%%%%%%%

             Acpt(iN, 1:size(cpt,1),1:size(cpt,2)) = cpt;
             Acpth(iN, 1:size(cpth,1),1:size(cpth,2)) = cpth;
             Awnt(iN, 1:size(wnt,1),1:size(wnt,2)) = wnt;
             Awnth(iN, 1:size(wnth,1),1:size(wnth,2)) = wnth;
             Acnt(iN, 1:size(cnt,1),1:size(cnt,2)) = cnt;
             Acnth(iN, 1:size(cnth,1),1:size(cnth,2)) = cnth;
             Awpt(iN, 1:size(wpt,1),1:size(wpt,2)) = wpt;
             Awpth(iN, 1:size(wpth,1),1:size(wpth,2)) = wpth;

             [Aa(iN, 1), b, c, d] = ttest2(mean(cnt(:,1318:1368),2), mean(wpt(:,1318:1368),2));
             [Aa(iN, 2), b, c, d] = ttest2(mean(cnth(:,1318:1368),2), mean(wpth(:,1318:1368),2));
             [Aa(iN, 3), b, c, d] = ttest2(mean(cpt(:,1318:1368),2), mean(wnt(:,1318:1368),2));
             [Aa(iN, 4), b, c, d] = ttest2(mean(cpth(:,1318:1368),2), mean(wnth(:,1318:1368),2));

        % % %      figure(940), clf, hold on,
        % % %      plot(mean(cpt(:,1:1700)), 'r.');plot(mean(wnt(:,1:1700)), 'k.');
        % % %      plot(mean(cpt(:,1:1700)) - mean(wnt(:,1:1700)),'r', 'LineWidth',2);
        % % %      
        % % %      plot(mean(cnt(:,1:1700)), 'g.');plot(mean(wpt(:,1:1700)), 'k-.');
        % % %      plot(mean(cnt(:,1:1700)) - mean(wpt(:,1:1700)),'g', 'LineWidth',2);
        % % %      
        % % %      plot(mean(cpth(:,1:1700)), 'b.');plot(mean(wnth(:,1:1700)), 'k-.');
        % % %      plot(mean(cpth(:,1:1700)) - mean(wnth(:,1:1700)),'b--', 'LineWidth',2);
        % % %      
        % % %      plot(mean(cnth(:,1:1700)), 'm.');plot(mean(wpth(:,1:1700)), 'k-.');
        % % %      %errorbar(mean(cnth(:,1:1700)), std(cnth(:,1:1700))./sqrt(size(cnth,1)), 'm.');
        % % %      plot(mean(cnth(:,1:1700)) - mean(wpth(:,1:1700)),'m--', 'LineWidth',2);
        % % %      
        % % %      xlim([0 1450])
        % % %      
        % % %      figure(640), clf, hold on,
        % % %      plot(mean(cptSacc),'r'); errorbar(mean(cptSacc),std(cptSacc)/sqrt(size(cptSacc,1)),'r');
        % % %      plot(mean(wntSacc), 'r.');
        % % %      plot(mean(wntSacc), 'b');
        % % %      plot(mean(cntSacc), 'm'); errorbar(mean(cntSacc),std(cntSacc)/sqrt(size(cntSacc,1)),'m');
        % % %      plot(mean(wptSacc), 'k.');
        % % %      plot(mean(wptSacc), 'k');
        % % %      xlim([-10 300]);
        % % %      
        % % %      stimes = []; sscc = {}; 
        % % %      for si = 1:length(saccads), 
        % % %          if ~isempty(saccads{si}), 
        % % %              sscc{end+1} = saccads{si}; 
        % % %              stimes(end+1) = saccads{si}.start; 
        % % %          end, 
        % % %      end
        % % %      figure(540), clf, 
        % % %      hist(stimes,20);
        % % %      
        % % %      stimes = []; sscc = {}; 
        % % %      for si = 1:length(saccads), 
        % % %          if ~isempty(saccads{si}) && abs(signal(si))>0.009, 
        % % %              sscc{end+1} = saccads{si}; 
        % % %              stimes(end+1) = saccads{si}.start; 
        % % %          end, 
        % % %      end
        % % %      hold on, hist(stimes,20, 'r');
        % % %      stimes = []; sscc = {}; 
        % % %      for si = 1:length(saccads), 
        % % %          if ~isempty(saccads{si}) && abs(signal(si))==0, 
        % % %              sscc{end+1} = saccads{si}; 
        % % %              stimes(end+1) = saccads{si}.start; 
        % % %          end, 
        % % %      end
        % % %      hold on, hist(stimes,20, 'r');
        % % %      
        % % %      h = findobj(gca,'Type','patch'); set(h,'FaceColor','r','EdgeColor','w')

             g1v=[];g2v=[];g3v=[];g4v=[];g5v=[];g6v=[];g7v=[];g8v=[];g9v=[];g10v=[];g11v=[];g12v=[];g13v=[];g14v=[];g15v=[];g16v=[];
             g17v=[];g18v=[];g19v=[];g20v=[];g21v=[];g22v=[];g23v=[];g24v=[];g25v=[];g26v=[];g27v=[];g28v=[];g29v=[];g30v=[];
             g1h=[];g2h=[];g3h=[];g4h=[];g5h=[];g6h=[];g7h=[];g8h=[];g9h=[];g10h=[];g11h=[];g12h=[];g13h=[];g14h=[];g15h=[];g16h=[];
             g17h=[];g18h=[];g19h=[];g20h=[];g21h=[];g22h=[];g23h=[];g24h=[];g25h=[];g26h=[];g27h=[];g28h=[];g29h=[];g30h=[];
             SacSpace = [];SColors = [];SSizes = [];
             for spi = 1:length(SacPool)
                 %SacSpace(spi, :, :, :) = [SacPool{spi}.start/10, SacPool{spi}.size, signal(SacPool{spi}.TrialIdx)];
                 %SacSpace(spi, :, :, :) = [SacPool{spi}.peakt/10, SacPool{spi}.size , signal(SacPool{spi}.TrialIdx)];
                 %SacSpace(spi, :, :, :) = [SacPool{spi}.peakt/10, SacPool{spi}.peakv , signal(SacPool{spi}.TrialIdx)];
                 %SacSpace(spi, :, :, :) = [SacPool{spi}.end/10, SacPool{spi}.peakv , signal(SacPool{spi}.TrialIdx)];
                 if SacPool{spi}.start > 20000 && SacPool{spi}.start < 25000 
                     %SacSpace(spi, :, :, :) = [(SacPool{spi}.end-SacPool{spi}.start)/10, SacPool{spi}.peakv , signal(SacPool{spi}.TrialIdx)];
                     %SacSpace(spi, :, :, :) = [(SacPool{spi}.end-SacPool{spi}.start)/10, abs(SacPool{spi}.dir) , signal(SacPool{spi}.TrialIdx)];
                     %SacSpace(spi, :, :, :) = [(SacPool{spi}.end-SacPool{spi}.start)/10, SacPool{spi}.size , signal(SacPool{spi}.TrialIdx)];
                     SacSpace(spi, :, :, :) = [SacPool{spi}.peakv, (SacPool{spi}.end-SacPool{spi}.start)/10 , signal(SacPool{spi}.TrialIdx)];
                     %disp(SacPool{spi}.TrialIdx);
                 else
                     SacSpace(spi, :, :, :) = [0 0 0];
                 end
                 sigsig = signal(SacPool{spi}.TrialIdx);
                 SColors(spi,:,:,:) = [0 0 1];
                 SSizes(spi) = 3;
                 %if sigsig == 0, 
                 %    SSizes(spi) = 40; SColors(spi,:,:,:) = [0 0 0];
                 %end
                 %if abs(sigsig) == max(signal), SColors(spi,:,:,:) = [1 0 0]; SSizes(spi) = 20; end


                 if cip0(SacPool{spi}.TrialIdx)  || cin0(SacPool{spi}.TrialIdx)  
                     SColors(spi,:,:,:) = [1 0.25 0]; SSizes(spi) = 43;
                     g1v(size(g1v,1)+1,:) = SacPool{spi}.lv;
                     g1h(size(g1h,1)+1,:) = SacPool{spi}.lh;
                 end
                 if win0(SacPool{spi}.TrialIdx)  || wip0(SacPool{spi}.TrialIdx)  
                     SColors(spi,:,:,:) = [1 0 1]; SSizes(spi) = 73;
                     g2v(size(g2v,1)+1,:) = SacPool{spi}.lv;
                     g2h(size(g2h,1)+1,:) = SacPool{spi}.lh;
                 end
                 if ~(cin0(SacPool{spi}.TrialIdx)  || wip0(SacPool{spi}.TrialIdx) || cip0(SacPool{spi}.TrialIdx)  || win0(SacPool{spi}.TrialIdx)  )
                     SSizes(spi) = 3;
                     g5v(size(g5v,1)+1,:) = SacPool{spi}.lv;
                     g5h(size(g5h,1)+1,:) = SacPool{spi}.lh;
                 end

                 if cip(SacPool{spi}.TrialIdx)  || cin(SacPool{spi}.TrialIdx)
                     SColors(spi,:,:,:) = [0 0.25 1]; SSizes(spi) = 23;
                     g3v(size(g3v,1)+1,:) = SacPool{spi}.lv;
                     g3h(size(g3h,1)+1,:) = SacPool{spi}.lh;
                 end
                 if win(SacPool{spi}.TrialIdx)  || wip(SacPool{spi}.TrialIdx)  
                     SColors(spi,:,:,:) = [0 1 1]; SSizes(spi) = 33;
                     g4v(size(g4v,1)+1,:) = SacPool{spi}.lv;
                     g4h(size(g4h,1)+1,:) = SacPool{spi}.lh;
                 end
                 %if ~(cin(SacPool{spi}.TrialIdx)  || wip(SacPool{spi}.TrialIdx) || cip(SacPool{spi}.TrialIdx)  || win(SacPool{spi}.TrialIdx)  )
                 %    SSizes(spi) = 3;
                 %end

                 if cip(SacPool{spi}.TrialIdx),  g11v(size(g11v,1)+1,:) = SacPool{spi}.lv; g11h(size(g11h,1)+1,:) = SacPool{spi}.lh; end
                 if cin(SacPool{spi}.TrialIdx),  g12v(size(g12v,1)+1,:) = SacPool{spi}.lv; g12h(size(g12h,1)+1,:) = SacPool{spi}.lh; end
                 if cip0(SacPool{spi}.TrialIdx), g13v(size(g13v,1)+1,:) = SacPool{spi}.lv; g13h(size(g13h,1)+1,:) = SacPool{spi}.lh; end
                 if cin0(SacPool{spi}.TrialIdx), g14v(size(g14v,1)+1,:) = SacPool{spi}.lv; g14h(size(g14h,1)+1,:) = SacPool{spi}.lh; end
                 if wip(SacPool{spi}.TrialIdx),  g15v(size(g15v,1)+1,:) = SacPool{spi}.lv; g15h(size(g15h,1)+1,:) = SacPool{spi}.lh; end
                 if win(SacPool{spi}.TrialIdx),  g16v(size(g16v,1)+1,:) = SacPool{spi}.lv; g16h(size(g16h,1)+1,:) = SacPool{spi}.lh; end
                 if wip0(SacPool{spi}.TrialIdx), g17v(size(g17v,1)+1,:) = SacPool{spi}.lv; g17h(size(g17h,1)+1,:) = SacPool{spi}.lh; end
                 if win0(SacPool{spi}.TrialIdx), g18v(size(g18v,1)+1,:) = SacPool{spi}.lv; g18h(size(g18h,1)+1,:) = SacPool{spi}.lh; end

             end

             if plotFigures
                 figure(333), clf, hold on,
                 if ~isempty(SacSpace)
                     scatter3(SacSpace(:,1), SacSpace(:,2), SacSpace(:,3), SSizes, SColors,'filled');
                     warning off
                     %xlim([1900 2600]);
                     %ylim([0.5 5]);
                     warning on
                     %zlim([]);
                 end     
                 figure(343), clf, hold on
                 plot(mean(g1v), 'r');
                 plot(mean(g2v), 'g');
                 plot(mean(g3v), 'b');
                 plot(mean(g4v), 'c');
                 plot(mean(g5v), 'm');

                 plot(mean(abs(g1v)), 'r.');
                 plot(mean(abs(g2v)), 'g.');
                 plot(mean(abs(g3v)), 'b.');
                 plot(mean(abs(g4v)), 'c.');
                 plot(mean(abs(g5v)), 'm.');
                 xlim([0 1667]);

                 figure(345), clf, hold on
                 plot(mean(abs(g11v)), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean(abs(g12v)), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean(abs(g13v)), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean(abs(g14v)), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean(abs(g15v)), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean(abs(g16v)), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean(abs(g17v)), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean(abs(g18v)), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1667]);
                 figure(346), clf, hold on
                 plot(mean(abs(g11h)), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean(abs(g12h)), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean(abs(g13h)), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean(abs(g14h)), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean(abs(g15h)), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean(abs(g16h)), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean(abs(g17h)), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean(abs(g18h)), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1667]);

                 figure(347), clf, hold on
                 plot(mean((g11v)), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean((g12v)), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean((g13v)), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean((g14v)), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean((g15v)), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean((g16v)), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean((g17v)), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean((g18v)), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1667]);
                 figure(348), clf, hold on
                 plot(mean((g11h)), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean((g12h)), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean((g13h)), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean((g14h)), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean((g15h)), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean((g16h)), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean((g17h)), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean((g18h)), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1667]);

                 figure(349), clf, hold on
                 plot(mean(cumsum(g11v')'), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean(cumsum(g12v')'), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean(cumsum(g13v')'), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean(cumsum(g14v')'), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean(cumsum(g15v')'), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean(cumsum(g16v')'), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean(cumsum(g17v')'), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean(cumsum(g18v')'), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1333]);%xlim([0 1667]);
                 figure(350), clf, hold on
                 plot(mean(cumsum(g11h')'), 'r.');  %plot(mean(g11v), 'r.');
                 plot(mean(cumsum(g12h')'), 'g.');  %plot(mean(g12v), 'g.');
                 plot(mean(cumsum(g13h')'), 'b.');  %plot(mean(g13v), 'b.');
                 plot(mean(cumsum(g14h')'), 'c.');  %plot(mean(g14v), 'c.');
                 plot(mean(cumsum(g15h')'), 'm.');  %plot(mean(g15v), 'm.');
                 plot(mean(cumsum(g16h')'), 'k.'); %plot(mean(g16v), 'k.-');
                 plot(mean(cumsum(g17h')'), 'k--');  %plot(mean(g17v), 'k.');
                 plot(mean(cumsum(g18h')'), 'r--'); %plot(mean(g18v), 'r.-');
                 xlim([0 1333]);%xlim([0 1667]);


                 figure(110), clf, hold on,
                 if ~isempty(g11v) && ~isempty(g17v) && ~isempty(g18v) && ~isempty(g15v)
                     scatter(mean(g11v(:,1:1333)), mean(g11h(:,1:1333)), 40, 'r');%, 'filled');
                     scatter(mean(g12v(:,1:1333)), mean(g12h(:,1:1333)), 40, 'g');%, 'filled');
                     scatter(mean(g13v(:,1:1333)), mean(g13h(:,1:1333)), 40, 'b');%, 'filled');
                     scatter(mean(g14v(:,1:1333)), mean(g14h(:,1:1333)), 40, 'c');%, 'filled');
                     scatter(mean(g15v(:,1:1333)), mean(g15h(:,1:1333)), 7, 'm', 'filled');
                     scatter(mean(g16v(:,1:1333)), mean(g16h(:,1:1333)), 7, 'y', 'filled');
                     scatter(mean(g17v(:,1:1333)), mean(g17h(:,1:1333)), 7, 'k', 'filled');
                     scatter(mean(g18v(:,1:1333)), mean(g18h(:,1:1333)), 7, [1 0.75 0.25], 'filled');
                     %xlim([-2.5 2.5]); ylim([-2.5 2.5])
                     xlim([-1.5 1.5]); ylim([-1.5 1.5])     %xlim([-.5 .5]); ylim([-.5 .5])
                 end
             end

             for di = 1:size(g11v,1), d11s(di) = ali2Ddistance([g11v(di,33:1000);g11h(di,33:1000)]); end, d11 = mean(d11s); d11se = std(d11s)/sqrt(length(d11s));   
             for di = 1:size(g12v,1), d12s(di) = ali2Ddistance([g12v(di,33:1000);g12h(di,33:1000)]); end, d12 = mean(d12s); d12se = std(d12s)/sqrt(length(d12s));  
             for di = 1:size(g13v,1), d13s(di) = ali2Ddistance([g13v(di,33:1000);g13h(di,33:1000)]); end, d13 = mean(d13s); d13se = std(d13s)/sqrt(length(d13s));  
             for di = 1:size(g14v,1), d14s(di) = ali2Ddistance([g14v(di,33:1000);g14h(di,33:1000)]); end, d14 = mean(d14s); d14se = std(d14s)/sqrt(length(d14s));  
             for di = 1:size(g15v,1), d15s(di) = ali2Ddistance([g15v(di,33:1000);g15h(di,33:1000)]); end, d15 = mean(d15s); d15se = std(d15s)/sqrt(length(d15s));  
             for di = 1:size(g16v,1), d16s(di) = ali2Ddistance([g16v(di,33:1000);g16h(di,33:1000)]); end, d16 = mean(d16s); d16se = std(d16s)/sqrt(length(d16s));  
             for di = 1:size(g17v,1), d17s(di) = ali2Ddistance([g17v(di,33:1000);g17h(di,33:1000)]); end, d17 = mean(d17s); d17se = std(d17s)/sqrt(length(d17s));  
             for di = 1:size(g18v,1), d18s(di) = ali2Ddistance([g18v(di,33:1000);g18h(di,33:1000)]); end, d18 = mean(d18s); d18se = std(d18s)/sqrt(length(d18s));  

             if plotFigures
                 figure(130), clf, hold on,
                 bar([d11, d12, d13, d14, d15, d16, d17, d18]);
                errorbar([d11, d12, d13, d14, d15, d16, d17, d18], [d11se, d12se, d13se, d14se, d15se, d16se, d17se, d18se]);
             end

             for deci = 1:size(g11v,1), dec11s(deci) = aliGazeEccentricity([g11v(deci,33:1000);g11h(deci,33:1000)]); end, dec11 = mean(dec11s); dec11se = std(dec11s)/sqrt(length(dec11s));   
             for deci = 1:size(g12v,1), dec12s(deci) = aliGazeEccentricity([g12v(deci,33:1000);g12h(deci,33:1000)]); end, dec12 = mean(dec12s); dec12se = std(dec12s)/sqrt(length(dec12s));  
             for deci = 1:size(g13v,1), dec13s(deci) = aliGazeEccentricity([g13v(deci,33:1000);g13h(deci,33:1000)]); end, dec13 = mean(dec13s); dec13se = std(dec13s)/sqrt(length(dec13s));  
             for deci = 1:size(g14v,1), dec14s(deci) = aliGazeEccentricity([g14v(deci,33:1000);g14h(deci,33:1000)]); end, dec14 = mean(dec14s); dec14se = std(dec14s)/sqrt(length(dec14s));  
             for deci = 1:size(g15v,1), dec15s(deci) = aliGazeEccentricity([g15v(deci,33:1000);g15h(deci,33:1000)]); end, dec15 = mean(dec15s); dec15se = std(dec15s)/sqrt(length(dec15s));  
             for deci = 1:size(g16v,1), dec16s(deci) = aliGazeEccentricity([g16v(deci,33:1000);g16h(deci,33:1000)]); end, dec16 = mean(dec16s); dec16se = std(dec16s)/sqrt(length(dec16s));  
             for deci = 1:size(g17v,1), dec17s(deci) = aliGazeEccentricity([g17v(deci,33:1000);g17h(deci,33:1000)]); end, dec17 = mean(dec17s); dec17se = std(dec17s)/sqrt(length(dec17s));  
             for deci = 1:size(g18v,1), dec18s(deci) = aliGazeEccentricity([g18v(deci,33:1000);g18h(deci,33:1000)]); end, dec18 = mean(dec18s); dec18se = std(dec18s)/sqrt(length(dec18s));  

             if plotFigures
                 figure(131), clf, hold on,
                 bar([d11, d12, d13, d14, d15, d16, d17, d18]);
                 errorbar([d11, d12, d13, d14, d15, d16, d17, d18], [d11se, d12se, d13se, d14se, d15se, d16se, d17se, d18se]);
             end 

             if ~isempty(g11v) && ~isempty(g17v) && ~isempty(g18v) && ~isempty(g15v)
                 d11vF = diff(g11v); d11vFm = mean(d11vF(:,1:1300)); d11vFse = std(d11vF(:,1:1300))/sqrt(length(d11vF(:,1:1300)));   
                 d12vF = diff(g12v); d12vFm = mean(d12vF(:,1:1300)); d12vFse = std(d12vF(:,1:1300))/sqrt(length(d12vF(:,1:1300)));   
                 d13vF = diff(g13v); d13vFm = mean(d13vF(:,1:1300)); d13vFse = std(d13vF(:,1:1300))/sqrt(length(d13vF(:,1:1300)));   
                 d14vF = diff(g14v); d14vFm = mean(d14vF(:,1:1300)); d14vFse = std(d14vF(:,1:1300))/sqrt(length(d14vF(:,1:1300)));   
                 d15vF = diff(g15v); d15vFm = mean(d15vF(:,1:1300)); d15vFse = std(d15vF(:,1:1300))/sqrt(length(d15vF(:,1:1300)));   
                 d16vF = diff(g16v); d16vFm = mean(d16vF(:,1:1300)); d16vFse = std(d16vF(:,1:1300))/sqrt(length(d16vF(:,1:1300)));   
                 d17vF = diff(g17v); d17vFm = mean(d17vF(:,1:1300)); d17vFse = std(d17vF(:,1:1300))/sqrt(length(d17vF(:,1:1300)));   
                 d18vF = diff(g18v); d18vFm = mean(d18vF(:,1:1300)); d18vFse = std(d18vF(:,1:1300))/sqrt(length(d18vF(:,1:1300)));   

                 d11hF = diff(g11h); d11hFm = mean(d11hF(:,1:1300)); d11hFse = std(d11hF(:,1:1300))/sqrt(length(d11hF(:,1:1300)));   
                 d12hF = diff(g12h); d12hFm = mean(d12hF(:,1:1300)); d12hFse = std(d12hF(:,1:1300))/sqrt(length(d12hF(:,1:1300)));   
                 d13hF = diff(g13h); d13hFm = mean(d13hF(:,1:1300)); d13hFse = std(d13hF(:,1:1300))/sqrt(length(d13hF(:,1:1300)));   
                 d14hF = diff(g14h); d14hFm = mean(d14hF(:,1:1300)); d14hFse = std(d14hF(:,1:1300))/sqrt(length(d14hF(:,1:1300)));   
                 d15hF = diff(g15h); d15hFm = mean(d15hF(:,1:1300)); d15hFse = std(d15hF(:,1:1300))/sqrt(length(d15hF(:,1:1300)));   
                 d16hF = diff(g16h); d16hFm = mean(d16hF(:,1:1300)); d16hFse = std(d16hF(:,1:1300))/sqrt(length(d16hF(:,1:1300)));   
                 d17hF = diff(g17h); d17hFm = mean(d17hF(:,1:1300)); d17hFse = std(d17hF(:,1:1300))/sqrt(length(d17hF(:,1:1300)));   
                 d18hF = diff(g18h); d18hFm = mean(d18hF(:,1:1300)); d18hFse = std(d18hF(:,1:1300))/sqrt(length(d18hF(:,1:1300)));   


                 d11vFmCS = cumsum(abs(d11vF)')'; d12vFmCS = cumsum(abs(d12vF)')'; d13vFmCS = cumsum(abs(d13vF)')';d14vFmCS = cumsum(abs(d14vF)')';
                 d15vFmCS = cumsum(abs(d15vF)')'; d16vFmCS = cumsum(abs(d16vF)')'; d17vFmCS = cumsum(abs(d17vF)')';d18vFmCS = cumsum(abs(d18vF)')';

                 d11vFmCSA = cumsum(cumsum(abs(d11vF)'))'; d12vFmCSA = cumsum(cumsum(abs(d12vF)'))'; d13vFmCSA = cumsum(cumsum(abs(d13vF)'))'; d14vFmCSA = cumsum(cumsum(abs(d14vF)'))';
                 d15vFmCSA = cumsum(cumsum(abs(d15vF)'))'; d16vFmCSA = cumsum(cumsum(abs(d16vF)'))'; d17vFmCSA = cumsum(cumsum(abs(d17vF)'))'; d18vFmCSA = cumsum(cumsum(abs(d18vF)'))';

                 d11vFmCSAa = diff(cumsum(abs(d11vF)')'); d12vFmCSAa = diff(cumsum(abs(d12vF)')'); d13vFmCSAa = diff(cumsum(abs(d13vF)')'); d14vFmCSAa = diff(cumsum(abs(d14vF)')');
                 d15vFmCSAa = diff(cumsum(abs(d15vF)')'); d16vFmCSAa = diff(cumsum(abs(d16vF)')'); d17vFmCSAa = diff(cumsum(abs(d17vF)')'); d18vFmCSAa = diff(cumsum(abs(d18vF)')');

                 d11hFmCSAa = diff(cumsum(abs(d11hF)')'); d12hFmCSAa = diff(cumsum(abs(d12hF)')'); d13hFmCSAa = diff(cumsum(abs(d13hF)')'); d14hFmCSAa = diff(cumsum(abs(d14hF)')');
                 d15hFmCSAa = diff(cumsum(abs(d15hF)')'); d16hFmCSAa = diff(cumsum(abs(d16hF)')'); d17hFmCSAa = diff(cumsum(abs(d17hF)')'); d18hFmCSAa = diff(cumsum(abs(d18hF)')');

             end

             if plotFigures
                 figure(141), clf, hold on
                 plot(d11vFm, 'r.');
                 plot(d12vFm, 'g.');
                 plot(d13vFm, 'b.');
                 plot(d14vFm, 'c.');
                 plot(d15vFm, 'm.');
                 plot(d16vFm, 'k.');
                 plot(d17vFm, 'k--');
                 plot(d18vFm, 'r--');


                 figure(142), clf, hold on,
                 bar([mean(d11vFm), mean(d12vFm), mean(d13vFm), mean(d14vFm), mean(d15vFm), mean(d16vFm), mean(d17vFm), mean(d18vFm)]);
                 %errorbar([mean(d11vFm), mean(d12vFm), mean(d13vFm), mean(d14vFm), mean(d15vFm), mean(d16vFm), mean(d17vFm), mean(d18)], [d11se, d12se, d13se, d14se, d15se, d16se, d17se, d18se]);

                 figure(143), clf, hold on
                 plot(mean(d11vFmCS), 'r.');
                 plot(mean(d12vFmCS), 'g.');
                 plot(mean(d13vFmCS), 'b.');
                 plot(mean(d14vFmCS), 'c.');
                 plot(mean(d15vFmCS), 'm.');
                 plot(mean(d16vFmCS), 'k.');
                 plot(mean(d17vFmCS), 'k--');
                 plot(mean(d18vFmCS), 'r--');
                 xlim([1 1556]);


                 figure(144), clf, hold on
                 plot(mean(d11vFmCSA), 'r.');
                 plot(mean(d12vFmCSA), 'g.');
                 plot(mean(d13vFmCSA), 'b.');
                 plot(mean(d14vFmCSA), 'c.');
                 plot(mean(d15vFmCSA), 'm.');
                 plot(mean(d16vFmCSA), 'k.');
                 plot(mean(d17vFmCSA), 'k--');
                 plot(mean(d18vFmCSA), 'r--');
                 xlim([1 1556]);

                 figure(145), clf, hold on
                 plot(mean(d11vFmCSAa), 'r.');
                 plot(mean(d12vFmCSAa), 'g.');
                 plot(mean(d13vFmCSAa), 'b.');
                 plot(mean(d14vFmCSAa), 'c.');
                 plot(mean(d15vFmCSAa), 'm.');
                 plot(mean(d16vFmCSAa), 'k.');
                 plot(mean(d17vFmCSAa), 'k--');
                 plot(mean(d18vFmCSAa), 'r--');
                 xlim([1 1556]);

                 figure(146), clf, hold on
                 plot(mean(d11hFmCSAa), 'r.');
                 plot(mean(d12hFmCSAa), 'g.');
                 plot(mean(d13hFmCSAa), 'b.');
                 plot(mean(d14hFmCSAa), 'c.');
                 plot(mean(d15hFmCSAa), 'm.');
                 plot(mean(d16hFmCSAa), 'k.');
                 plot(mean(d17hFmCSAa), 'k--');
                 plot(mean(d18hFmCSAa), 'r--');
                 xlim([1 1556]);   
             end     

             figure(147), clf, hold on
             plot(abs(mean(d11hFmCSAa))+abs(mean(d11hFmCSAa)), 'r.');
             plot(abs(mean(d12hFmCSAa))+abs(mean(d12hFmCSAa)), 'g.');
             plot(abs(mean(d13hFmCSAa))+abs(mean(d13hFmCSAa)), 'b.');
             plot(abs(mean(d14hFmCSAa))+abs(mean(d14hFmCSAa)), 'c.');
             plot(abs(mean(d15hFmCSAa))+abs(mean(d15hFmCSAa)), 'm.');
             plot(abs(mean(d16hFmCSAa))+abs(mean(d16hFmCSAa)), 'k.');
             plot(abs(mean(d17hFmCSAa))+abs(mean(d17hFmCSAa)), 'k--');
             plot(abs(mean(d18hFmCSAa))+abs(mean(d18hFmCSAa)), 'r--');
             xlim([1 1556]);   
             title(num2str(iN));
             print(gcf, '-dtiff', strcat('/Users/ali/Desktop/figs/f147-',num2str(iN),'.tiff'));

             dXXhFmCSAa{iN} = {d11hFmCSAa, d12hFmCSAa, d13hFmCSAa, d14hFmCSAa, d15hFmCSAa, d16hFmCSAa, d17hFmCSAa, d18hFmCSAa};

             if ~isempty(g11v)
                 figure(777), clf, hold on,
                 try
                 v11v = mean(var(g11v(:,1:1334)'));v12v = mean(var(g12v(:,1:1334)'));v13v = mean(var(g13v(:,1:1334)'));v14v = mean(var(g14v(:,1:1334)'));
                 v15v = mean(var(g15v(:,1:1334)'));v16v = mean(var(g16v(:,1:1334)'));v17v = mean(var(g17v(:,1:1334)'));v18v = mean(var(g18v(:,1:1334)'));
                 v11h = mean(var(g11h(:,1:1334)'));v12h = mean(var(g12h(:,1:1334)'));v13h = mean(var(g13h(:,1:1334)'));v14h = mean(var(g14h(:,1:1334)'));
                 v15h = mean(var(g15h(:,1:1334)'));v16h = mean(var(g16h(:,1:1334)'));v17h = mean(var(g17h(:,1:1334)'));v18h = mean(var(g18h(:,1:1334)'));

                 bar([v11v, v12v, v13v, v14v, v15v, v16v, v17v, v18v; v11h, v12h, v13h, v14h, v15h, v16h, v17h, v18h]'); 
                 set(gca,'XTickLabel',{'cip', 'cin', 'cip0', 'cin0', 'wip', 'win', 'wip0', 'win0'}');

                 figure(778), clf, hold on,
                 bar([v11v, v16v, v12v, v15v, v13v, v18v, v14v, v17v; v11h, v16h, v12h, v15h v13h, v18h, v14h, v17h]'); 
                 %set(gca,'XTickLabel',{'cip', 'cin', 'cip0', 'cin0', 'wip', 'win', 'wip0', 'win0'}');
                 catch exception
                     disp(exception.message);
                 end
             end

             SacMetr = [];
             for sci = 1:length(AllSacc)
                 if ~isempty(AllSacc{sci})
                     SacMetr(sci,:,:) = [length(AllSacc{sci}), mean([AllSacc{sci}([AllSacc{sci}.start]<20000).size]), sum([AllSacc{sci}([AllSacc{sci}.start]<20000).size])];
                 end
             end
             if plotFigures
                 figure(876), clf, hold on
                 scatter(SacMetr(cip,3), SacMetr(cip,2), 'r', 'filled')
                 scatter(SacMetr(win,3), SacMetr(win,2), 'b')
                 scatter(SacMetr(cin,3), SacMetr(cin,2), 'r', 'filled')
                 scatter(SacMetr(wip,3), SacMetr(wip,2), 'b')
                 scatter(SacMetr(cip0,3), SacMetr(cip0,2), 'g', 'filled')
                 scatter(SacMetr(win0,3), SacMetr(win0,2), 'c', 'filled')
                 scatter(SacMetr(cin0,3), SacMetr(cin0,2), 'm', 'filled')
                 scatter(SacMetr(wip0,3), SacMetr(wip0,2), 'k', 'filled')
             end     
             debugpoint = 1;
     end
end

%%

figure(888), clf, hold on
plot(squeeze(mean(mean(Acnt(9:39,:,:),2),1)) - squeeze(mean(mean(Awpt(9:39,:,:),2),1)),'r', 'LineWidth',2)
plot(squeeze(mean(mean(Acnth(9:39,:,:),2),1)) - squeeze(mean(mean(Awpth(9:39,:,:),2),1)),'g', 'LineWidth',2)
plot(squeeze(mean(mean(Acpt(9:39,:,:),2),1)) - squeeze(mean(mean(Awnt(9:39,:,:),2),1)),'b--', 'LineWidth',2)
plot(squeeze(mean(mean(Acpth(9:39,:,:),2),1)) - squeeze(mean(mean(Awnth(9:39,:,:),2),1)),'m--', 'LineWidth',2)
xlim([0 1500])

figure(891), clf, hold on
plot((squeeze(mean(mean(Acnt(9:39,:,1:1700),2),1)) - squeeze(mean(mean(Awpt(9:39,:,1:1700),2),1))) + (squeeze(mean(mean(Acnth(9:39,:,1:1700),2),1)) - squeeze(mean(mean(Awpth(9:39,:,1:1700),2),1))),'g', 'LineWidth',2);
plot((squeeze(mean(mean(Acpt(9:39,:,1:1700),2),1)) - squeeze(mean(mean(Awnt(9:39,:,1:1700),2),1))) + (squeeze(mean(mean(Acpth(9:39,:,1:1700),2),1)) - squeeze(mean(mean(Awnth(9:39,:,1:1700),2),1))),'m--', 'LineWidth',2)
xlim([0 1500])

figure(892), clf, hold on
plot((squeeze(mean(mean(Acnt(sum(Aa,2)>0,:,1:1700),2),1)) - squeeze(mean(mean(Awpt(sum(Aa,2)>0,:,1:1700),2),1))) + (squeeze(mean(mean(Acnth(sum(Aa,2)>0,:,1:1700),2),1)) - squeeze(mean(mean(Awpth(sum(Aa,2)>0,:,1:1700),2),1))),'g', 'LineWidth',2);
plot((squeeze(mean(mean(Acpt(sum(Aa,2)>0,:,1:1700),2),1)) - squeeze(mean(mean(Awnt(sum(Aa,2)>0,:,1:1700),2),1))) + (squeeze(mean(mean(Acpth(sum(Aa,2)>0,:,1:1700),2),1)) - squeeze(mean(mean(Awnth(sum(Aa,2)>0,:,1:1700),2),1))),'m--', 'LineWidth',2)
xlim([0 1500])

figure(893), clf, hold on
plot((squeeze(mean(mean(Acnt(sum(Aa,2)>1,:,1:1700),2),1)) - squeeze(mean(mean(Awpt(sum(Aa,2)>1,:,1:1700),2),1))) + (squeeze(mean(mean(Acnth(sum(Aa,2)>1,:,1:1700),2),1)) - squeeze(mean(mean(Awpth(sum(Aa,2)>1,:,1:1700),2),1))),'g', 'LineWidth',2);
plot((squeeze(mean(mean(Acpt(sum(Aa,2)>1,:,1:1700),2),1)) - squeeze(mean(mean(Awnt(sum(Aa,2)>1,:,1:1700),2),1))) + (squeeze(mean(mean(Acpth(sum(Aa,2)>1,:,1:1700),2),1)) - squeeze(mean(mean(Awnth(sum(Aa,2)>1,:,1:1700),2),1))),'m--', 'LineWidth',2)
xlim([0 1500])

figure(894), clf, hold on
plot((squeeze(mean(mean(Acnt(sum(Aa,2)>2,:,1:1700),2),1)) - squeeze(mean(mean(Awpt(sum(Aa,2)>2,:,1:1700),2),1))) + (squeeze(mean(mean(Acnth(sum(Aa,2)>2,:,1:1700),2),1)) - squeeze(mean(mean(Awpth(sum(Aa,2)>2,:,1:1700),2),1))),'g', 'LineWidth',2);
plot((squeeze(mean(mean(Acpt(sum(Aa,2)>2,:,1:1700),2),1)) - squeeze(mean(mean(Awnt(sum(Aa,2)>2,:,1:1700),2),1))) + (squeeze(mean(mean(Acpth(sum(Aa,2)>2,:,1:1700),2),1)) - squeeze(mean(mean(Awnth(sum(Aa,2)>2,:,1:1700),2),1))),'m--', 'LineWidth',2)
xlim([0 1500])


figure(895), clf, hold on
plot((squeeze(mean(mean(Acnt(sum(Aa,2)>3,:,1:1700),2),1)) - squeeze(mean(mean(Awpt(sum(Aa,2)>3,:,1:1700),2),1))) + (squeeze(mean(mean(Acnth(sum(Aa,2)>3,:,1:1700),2),1)) - squeeze(mean(mean(Awpth(sum(Aa,2)>3,:,1:1700),2),1))),'g', 'LineWidth',2);
plot((squeeze(mean(mean(Acpt(sum(Aa,2)>3,:,1:1700),2),1)) - squeeze(mean(mean(Awnt(sum(Aa,2)>3,:,1:1700),2),1))) + (squeeze(mean(mean(Acpth(sum(Aa,2)>3,:,1:1700),2),1)) - squeeze(mean(mean(Awnth(sum(Aa,2)>3,:,1:1700),2),1))),'m--', 'LineWidth',2)
xlim([0 1500])

