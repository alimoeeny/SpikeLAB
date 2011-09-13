function [rv] = expMining(Expt, Trials, Prop)

if isempty(Trials)
    Trials = 1:length([Expt.Trials]);
end

r = [];
for i = 1: size(Trials, 1)
    r(i, :) = [Expt.Trials(Trials(i,:)).(Prop)];
end

rv = r;
