validCells = []; for i = 1: length(PSTHs), if~isempty(PSTHs{i}), validCells(i) = [PSTHs{i}{4}]; end; end
AllConditions = [];

for i = 1 : length(AllNeurons),
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
        ebSacTrig = ([PSTHs{i}{13}]);
        for iST = 1: size(ebSacTrig,1),
            ebSacTrigMean(iST,:) = squeeze(mean(ebSacTrig(iST,squeeze(mean(ebSacTrig(iST,:,:),3))'>0,:)));
        end
        tPSTHSacTrigMean(i, :, :) = ebSacTrigMean;
                
%         figure(11), clf, hold on,
%         for tmpi = 1:30,
%             plot(ebSacTrigMean(tmpi,:)); 
%         end
        
        
        
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