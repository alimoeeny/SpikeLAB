load ~/Desktop/matlab.cylDPI.new.mat

m = 0;
for iN = 1:88
    try
        load(['~/Desktop/matlab.',num2str(iN),'.mat']);
        AllPSTHs{iN} = PSTHs{iN};
        AllTIs(iN) = TI(iN);
    catch
        disp(['Missing ', num2str(iN) ]);
        m = m + 1;
    end
end
disp(['Total Missing: ' , num2str(m), ' out of 88 that is: ', num2str(m/88)]);



%%
PSTHs  = AllPSTHs;
TI = AllTIs;

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
    close all;
    figure(2502),
    clf, hold on,
    a = [];
    h = plot(squeeze(mean(squeeze(tPSTHSacTrigMean(sum(tPSTHSacTrigMean(:, 41,:),3)~=0,41,:)))), 'r');
    set(h, 'LineWidth', 6);
    h = plot(squeeze(mean(squeeze(tPSTHSacTrigMean(sum(tPSTHSacTrigMean(:, 42,:),3)~=0,42,:)))), 'b');
    set(h, 'LineWidth', 6);
    set(gca, 'XGrid', 'on');
    set(gca, 'GridLineStyle','--');
    xlim([0 480]);
    xtl = [-100, 0, 50, 100, 200, 300];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    title(FileType);

    
    
    
    
%%
figure(8392), clf , hold on,
% a(:,1) = squeeze(mean(tPSTHSacTrigMean(TI>0,41,:)));
% a(:,2) = squeeze(mean(tPSTHSacTrigMean(TI<0,41,:)));
% 
% a(:,3) = squeeze(mean(tPSTHSacTrigMean(TI>0,42,:)));
% a(:,4) = squeeze(mean(tPSTHSacTrigMean(TI<0,42,:)));

a(:,1) = squeeze(mean(tPSTHSacTrigMean(PIS(1:size(tPSTHSacTrigMean,1),10)>0,41,:)));
a(:,2) = squeeze(mean(tPSTHSacTrigMean(PIS(1:size(tPSTHSacTrigMean,1),10)<0,41,:)));

a(:,3) = squeeze(mean(tPSTHSacTrigMean(PIS(1:size(tPSTHSacTrigMean,1),10)>0,42,:)));
a(:,4) = squeeze(mean(tPSTHSacTrigMean(PIS(1:size(tPSTHSacTrigMean,1),10)<0,42,:)));


h = plot(mean(a(:,[1]),2), 'r');
set(h, 'LineWidth',6);
h = plot(mean(a(:,[2]),2), 'r');
set(h, 'LineWidth',6);
h = plot(mean(a(:,[3]),2), 'b');
set(h, 'LineWidth',6);
h = plot(mean(a(:,[4]),2), 'b');
set(h, 'LineWidth',6);

h = plot(mean(a(:,[1,4]),2), 'r');
set(h, 'LineWidth',2);
h = plot(mean(a(:,[2,3]),2), 'b');
set(h, 'LineWidth',2);

set(gca, 'XGrid', 'on');
xlim([0 480]);
xtl = [-100, 0, 50, 100, 200, 300];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
set(gca, 'GridLineStyle','--');
title(FileType);
