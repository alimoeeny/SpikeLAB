
%%
load ~/Desktop/DistanceFileNames DistanceFileNames

for i = 1 : length(DistanceFileNames)
    dm{i} = load(DistanceFileNames{i});
end

%%

for i = 1:length(dm)
    DM = dm{i}.DistanceMatrix{1}.D(:,:,1);

    [scores, permutations] = forceSqueezeDistanceMatrix(DM);

    maxidx = find(scores==max(scores));
    figure(959), clf, hold on, 
    set(gcf, 'PaperPositionMode', 'auto')
    set(gcf, 'Position', [100 200 800 250]);
    
    subplot(1,3,1);
    surf(tril(DM([1:end end], [1:end end])));
    view(2);
    set(gca, 'CLim', [0 5]);
    set(gca, 'Xlim', [1 7])
    
    remapedDM = remapDistanceMatrix(DM, permutations(maxidx(1),:));
    subplot(1,3,2);
    surf(tril(remapedDM([1:end end], [1:end end])));
    view(2);
    set(gca, 'CLim', [0 5]);
    set(gca, 'Xlim', [1 7])
    
    if(length(maxidx)>1)
        remapedDM = remapDistanceMatrix(DM, permutations(maxidx(2),:));
        subplot(1,3,3);
        surf(tril(remapedDM([1:end end], [1:end end])));
        view(2);
        set(gca, 'CLim', [0 5]);
        set(gca, 'Xlim', [1 7])
    end
    print(959, '-depsc', '-tiff', '-r300', [num2str(i), 'distmatgrid.eps']);

end