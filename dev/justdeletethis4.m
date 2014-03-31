AvgWtPSTH = [];
for wc = 1: size(AllConditions,2)
    disp(wc);
    avgbuff = [];
    for wi =1: size(tPSTH,1)
        if (AllConditions(wi,wc)>0)
            avgbuff = [avgbuff, squeeze(tPSTH(wi,wc,:))];
        end
    end
    AvgWtPSTH(wc,:) = mean(avgbuff');
end
