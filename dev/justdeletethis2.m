for i = 1: length(AllNeurons),


    figure(2502),
    clf, hold on,
    a = [];
    a(:,1) = squeeze(tPSTHSacTrigMean(i,41,:));
    a(:,2) = squeeze(tPSTHSacTrigMean(i,42,:));
    h = plot(a);
    set(h, 'LineWidth', 2);
    set(gca, 'XGrid', 'on');
    xlim([0 500]);
    xtl = [-100, 0, 50, 100, 250, 500];
    set(gca, 'XTick', xtl+100-(BinSize - SmoothingBinSize)/2);
    set(gca, 'XTickLabel', {num2str(xtl')});
    legend(h, GetLegends(FileType));
    title(FileType);
end
