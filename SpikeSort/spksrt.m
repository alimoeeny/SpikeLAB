clear


%%
templates = [];
tmpT = load('spikeTemplate000.mat');
templates(1,:) = tmpT.rndspk;
tmpT = load('spikeTemplate001.mat');
templates(2,:) = tmpT.msnorm;
tmpT = load('spikeTemplate002.mat');
templates(3,:) = tmpT.msnorm;

%%
%load /sd/bgc6/smr/dae/M002/Expt6FullV.mat
load /data/Expt6FullV.mat

%%
p6 = cast(FullV.V(6,:), 'double');
p5 = FullV.V(5,:);
p7 = FullV.V(7,:);

p75 = (p7 + p5)./2;

%%
windowSize = 25;
threshold = 10000;
 
%%
[swatches, timestamps] = spikesOut(p6, threshold, windowSize);

%%

p6norm = (p6 - mean(p6)) ./ std(p6);

threshold = 4;
[swatchesnorm, timestampsnorm] = spikesOut(p6norm, threshold, windowSize);


%%
templateSpace = [];

for ti = 1 : size(templates,1)
    for i = 1: size(swatchesnorm,1)
        templateSpace(ti, i) = sum(templates(ti,:) * swatchesnorm(i));
    end
end

%%

ms = mean(swatches);

p6Cms = conv(cast(p6, 'double'), cast(ms, 'double'), 'full');
%p6Cms = (p6Cms - mean(p6Cms)) ./ range(p6Cms);
p6Cms = (p6Cms - mean(p6Cms)) ./ std(p6Cms);

%%
threshold = 4.2;
[mcswatches, mctimestamps] = spikesOut(p6Cms, threshold, windowSize);

figure(1878), scatter(sum(abs(mcswatches)), var(mcswatches))

%%
for i = 1: size(swatches,1)
   templatescores(i) = sum(abs(ms + swatches(i,:))); 
end

