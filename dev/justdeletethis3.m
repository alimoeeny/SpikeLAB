for iN = 1:88
    try
        load(['~/Desktop/matlab.',num2str(iN),'.mat']);
        AllPSTHs{iN} = PSTHs{iN};
        AllTIs(iN) = TI(iN);
    catch
        disp(['Missing ', num2str(iN) ]);
    end
end




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
    figure(2502),
    clf, hold on,
    a = [];
    h = plot(squeeze(mean(squeeze(tPSTHSacTrigMean(sum(tPSTHSacTrigMean(:, 41,:),3)~=0,41,:)))));
    set(h, 'LineWidth', 6);
    h = plot(squeeze(mean(squeeze(tPSTHSacTrigMean(sum(tPSTHSacTrigMean(:, 42,:),3)~=0,42,:)))), 'r');
    set(h, 'LineWidth', 6);
    set(gca, 'XGrid', 'on');
    xlim([0 480]);
    xtl = [-100, 0, 50, 100, 200, 300];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    title(FileType);
