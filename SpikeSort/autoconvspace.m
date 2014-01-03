load /data/Expt6FullV.mat
%%
V12 = cast(FullV.V(12,:), 'double') ;
V12n = V12 ./ (2 * std(V12));

%%
winsize = 30;

spansize = min(10000000, length(V12n));
convscore = [];
streamoffset = 1100000;
%par
parfor i = 1:100000 %length(V12n)-winsize
    step = i; % * winsize/3;
    sliwin = V12n(step+ streamoffset:step + streamoffset+winsize);
    convscore(i) = sum(conv(V12n(1:spansize), sliwin, 'same')./ sum(sliwin));
end

%disp('doing the normalization')
%convscoreN = convscore ./ (2 * std(convscore));



%%
disp('doing the spikesOut thing now');
[s, t] = spikesOutResearch(V12n(streamoffset:streamoffset+length(convscore)-1), convscore, 1, 30);
