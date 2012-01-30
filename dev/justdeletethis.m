for i = 1: length(AllNeurons),
    PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean(1:i,41:42,:)));
    figure(007007), clf, hold on,
    h = plot(squeeze(PopPSTHSacTrigMean)');
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
    xlim([0 500]);
    xtl = [-100, 0, 50, 100, 250, 500];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    legend(h, GetLegends(FileType));
    title(FileType);
end


%%
for i = 1: length(AllNeurons),
    PopPSTHSacTrigMean = squeeze((tPSTHSacTrigMean(i,[41:42],:)));
    figure(007007), clf, hold on,
    h = plot(squeeze(PopPSTHSacTrigMean)');
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
    xlim([0 500]);
    xtl = [-100, 0, 50, 100, 250, 500];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    legend(h, GetLegends(FileType));
    title(FileType);
end


%%
PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean(PIS(:,10)>0,41:end,:)));
figure(008008), clf, hold on,
h = plot(squeeze(PopPSTHSacTrigMean)');
PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean(PIS(:,10)<0,41:end,:)));
hold on, 
plot(squeeze(PopPSTHSacTrigMean)');
set(h, 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);


%%
selector = [1:1:34];
selector(88)=0;
PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean((TI>0) & (selector > 0) , 41:42,:)));
figure(010010), clf, hold on,
h = plot(squeeze(PopPSTHSacTrigMean)');
PopPSTHSacTrigMean = squeeze(mean(tPSTHSacTrigMean((TI<0) & (selector > 0) , 41:42,:)));
hold on, 
plot(squeeze(PopPSTHSacTrigMean)');
set(h, 'LineWidth', 2);
set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);


%%
figure(8392), clf , hold on,
a(:,1) = squeeze(mean(tPSTHSacTrigMean(TI>0,41,:)));
a(:,2) = squeeze(mean(tPSTHSacTrigMean(TI<0,41,:)));

a(:,3) = squeeze(mean(tPSTHSacTrigMean(TI>0,42,:)));
a(:,4) = squeeze(mean(tPSTHSacTrigMean(TI<0,42,:)));

h = plot(mean(a(:,[1]),2));
h = plot(mean(a(:,[2]),2));
h = plot(mean(a(:,[3]),2));
h = plot(mean(a(:,[4]),2));
h = plot(mean(a(:,[1,4]),2));
h = plot(mean(a(:,[2,3]),2));

set(h, 'LineWidth',2);
set(gca, 'XGrid', 'on');
xlim([0 500]);
xtl = [-100, 0, 50, 100, 250, 500];
set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
set(gca, 'XTickLabel', {num2str(xtl')});
legend(h, GetLegends(FileType));
title(FileType);
