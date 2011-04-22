function c = NotDoingVeryGood(Expt)

values = unique([Expt.Trials(:).dx]);
values = values(values<0.5 & values>-0.5);
tmpdxvals = [];
for dxx = 1: length(values)
   if (sum([Expt.Trials(:).dx]==values(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4 && sum([Expt.Trials(:).dx]== -values(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4)
       tmpdxvals(end+1) = values(dxx);
   end
end
values = tmpdxvals;     


if iseven(length(values))
    disp('WHAT IS GOING ON HERE?')
end

for i = 1:length(values)
    r(i) = sum([Expt.Trials(:).dx]==values(i) & [Expt.Trials(:).RespDir]==1)/sum([Expt.Trials(:).dx]==values(i) & [Expt.Trials(:).RespDir]~=0);
end

for i = 1:length(Expt.Trials)
    if sum(values==[Expt.Trials(i).dx])==0
        c(i) = 0;
    else
        c(i) = r(values==[Expt.Trials(i).dx])>=0.20 & r(values==[Expt.Trials(i).dx])<=0.80;
    end
end


