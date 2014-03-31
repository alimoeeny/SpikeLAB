% streamliner

clear
clc

tmp = load('../AllFullVFiles.mat', 'AllFullVFiles');
AllFullVFiles  = tmp.AllFullVFiles;
clear tmp;

gaussCount = 9;


%par
for iFV = 1: length(AllFullVFiles) % 25101
    slashes = strfind(AllFullVFiles{iFV}, '/');
    filePath = ['/data' AllFullVFiles{iFV}(slashes(4):end)];
    disp(filePath);
    tmp = load(filePath);
    FullV = tmp.FullV;
    clear tmp;
    
     figure(2243+iFV), clf, hold on
     set(gcf, 'Position', [2500 700 1200 1200]);
     figure(1166+iFV), clf, hold on
     set(gcf, 'Position', [1450 700 1200 1200]);
     figure(1819+iFV), clf, hold on
     set(gcf, 'Position', [300 700 1200 1200]);
     figure(1212+iFV), clf, hold on
     set(gcf, 'Position', [1 700 1200 1200]);
     
    for vi = 1:size(FullV.V,1)
        [timestamps] = detectEventsinRawVoltage(FullV.V, vi, 50, 'down');
        swatches = chopRawVoltage(FullV.V, vi, timestamps);
        %[coeff,score,latent,tsquared,explained] = pca(swatches);
        
        dSWtmp = cell(size(swatches,1),1);
        %par
        parfor iSW = 1:size(swatches,1)
            dswtemp = zeros(size(swatches,1),1);
            for jSW = iSW:size(swatches,1)
                %dswtemp(jSW,:) = sum(conv(swatches(iSW,:), swatches(jSW, :), 'same')) ./ (sum(swatches(iSW,:)) + sum(swatches(jSW, :)));
                [dswtemp(jSW,:), ~] = corr(swatches(iSW,:)', swatches(jSW, :)');
            end
            dSWtmp{iSW} = sign(dswtemp) .* (1 - abs(dswtemp));
        end
        dSW = zeros(size(swatches,1), size(swatches,1));
        for iSW = 1:size(swatches,1)
            dSW(:, iSW) = dSWtmp{iSW}';
        end
        
        %dSW = dSW + flipud(rot90(dSW));
        for junki = 1:size(dSW,1), dSW(junki,junki)=0; end        
        
        sqd = squareform(dSW(1:100,1:100))';
        lk = linkage(sqd, 'single');
        T = cluster(lk, 'maxclust', 10);

        
        spikeTemplates = LoadAllSpikeTemplates(30);
        TScoresCell = {};
        parfor iST = 1:size(spikeTemplates,1)
            tempscores = [];
            for iSW = 1:size(swatches,1)
                 tempscores(iSW,:) = sum(conv(spikeTemplates(iST,:), swatches(iSW,:), 'same')) ./ (2 * std(swatches(iSW,:)));
            end
            TScoresCell{iST} = tempscores;
        end
        TScores = zeros(size(spikeTemplates,1),1);
        for iST = 1:size(spikeTemplates,1)
            %for iSW = 1:size(swatches,1)
                TScores(iST, :) = TScoresCell{iST};
            %end
        end
        
        %AllScores = [TScores' score];
        [coeff,score,latent,tsquared,explained] = pca(TScores');
        AllScores = score;
        
        figure(1819+iFV),
        subplot(6,4,vi); hold on,
        %X = score(:,1:4);
        X = AllScores(:,1:4); %][1 7 9 10 11 12 13 14 15 16]);
        %scatter(X(:,1),X(:,2),10,'ko')
        options = statset('Display','final', 'UseParallel', 'always', 'MaxIter', 90);
        gm = gmdistribution.fit(X,gaussCount,'Options',options);
        %ezcontour(@(x,y)pdf(gm,[x y]),[-8 6],[-8 6]);
        idx = cluster(gm,X);
        %cluster1 = (idx == 1);
        %cluster2 = (idx == 2);
        %cluster3 = (idx == 3);
        %scatter(X(cluster1,1),X(cluster1,2),10,'r+');
        %scatter(X(cluster2,1),X(cluster2,2),10,'bo');
        for iIdx = 1:gaussCount
            ci = mod([iIdx / 18 iIdx / 6 iIdx / gaussCount], 1);
            scatter(X((idx == iIdx),1),X((idx==iIdx),2),10,ci);
        end
        
        
        [idxk,C,sumd,D] = kmeans(X,gaussCount,'Distance','city',...
                    'Replicates',5,...
                    'Options',options);
        
        figure(1212+iFV),
        subplot(6,4,vi); hold on,
        for iIdx = 1:gaussCount
            ci = mod([iIdx / 18 iIdx / 6 iIdx / gaussCount], 1);
            scatter(X((idxk == iIdx),1),X((idxk==iIdx),2),10,ci);
        end
        
         DMk = zeros(gaussCount, gaussCount);
         for iIdx = 1:gaussCount
             for jIdx = 1:gaussCount
                 DMk(iIdx, jIdx) = sqrt(sum((C(iIdx,:) - C(jIdx,:)).^2)); %mean(MahalDs(idx==iIdx, jIdx));
                 %DMk(iIdx, jIdx) = DMk(iIdx, jIdx) ./ (sqrt(sum(sum(gm.Sigma(:,:,iIdx))) + sum(sum(gm.Sigma(:,:,jIdx))))./1.4142);
             end
         end
         
        [squishedMatrix, score, cellMap] = squishDistanceMatrix( DMk );
        figure(2243+iFV), 
        subplot(6,4, vi); hold on
        %imagesc(squishedMatrix);
        sqdm = tril(squishedMatrix([1:end end], [1:end end]));
        surf(sqdm);
        view(2);
        colorbar('peer', gca);
        %set(gca, 'CLim', [0 5]);
        set(gca, 'Xlim', [1 gaussCount+1])
        set(gca, 'Ylim', [1 gaussCount+1])
        
        figure(1166+iFV),
        subplot(6,4,vi); hold on
        for iIdx = 1:gaussCount
            ci = mod([iIdx / 18 iIdx / 6 iIdx / gaussCount], 1);
            plot(swatches((idxk == iIdx),:)', 'Color', ci);
        end
        
        
        
        
%         MahalDs = mahal(gm, X); 
%         DM = [];
%         for iIdx = 1:gaussCount
%             for jIdx = 1:gaussCount
%                 DM(iIdx, jIdx) = sqrt(sum((gm.mu(iIdx,:) - gm.mu(jIdx,:)).^2)); %mean(MahalDs(idx==iIdx, jIdx));
%                 DM(iIdx, jIdx) = DM(iIdx, jIdx) ./ (sqrt(sum(sum(gm.Sigma(:,:,iIdx))) + sum(sum(gm.Sigma(:,:,jIdx))))./1.4142);
%             end
%         end
%         
%         [squishedMatrix, score, cellMap] = squishDistanceMatrix( DM );
%         figure(2243+iFV), 
%         subplot(6,4, vi); hold on
%         %imagesc(squishedMatrix);
%         sqdm = tril(squishedMatrix([1:end end], [1:end end]));
%         surf(sqdm);
%         view(2);
%         colorbar('peer',gca);
%         set(gca, 'CLim', [0 5]);
%         set(gca, 'Xlim', [1 gaussCount+1])
%         set(gca, 'Ylim', [1 gaussCount+1])
    end
    
    disp('PRESS A KEY WHEN READY TO DO THE NEXT ...');
    pause
    %if(ishandle(1819+iFV-1))
    %    close(1819+iFV-1);
    %end

end