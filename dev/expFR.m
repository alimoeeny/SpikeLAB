function [rv] = expFR(Expt, StartTime, FinishTime, varargin)


if length(varargin) < 2
    SpikeCounts = zeros(length([Expt.Trials]),1);
    for tr = 1: length([Expt.Trials]), 
        SpikeCounts(tr) = sum([Expt.Trials(tr).Spikes]>=StartTime & [Expt.Trials(tr).Spikes]<=FinishTime);
    end
    fr = SpikeCounts * (1000 / ((FinishTime - StartTime) / 10 ));
    
    rv = fr;
    return
else
    allFRs = expFR(Expt, StartTime, FinishTime);
    Props = [];
    Values = {};
    for i = 1:length(varargin)
        Props(i, :) = varargin{i};
        Values{i} = unique(Props(i,:));
    end

    conditionsSize = [];
    for i = 1:length(Values)
        conditionsSize = [conditionsSize, size(Values{i},2)];
    end

    FRs = zeros(conditionsSize);
    FSEs= zeros(conditionsSize);
    for i = 1:size(Values{1},2)
       for j = 1:size(Values{2},2)
           Trials = (Props(1,:) == Values{1}(i)) & (Props(2,:) == Values{2}(j));
           if(sum(Trials)==0)
               FRs(i,j) = -0.123456789;
               FSEs(i,j)= -0.123456789;
           else
               FRs(i,j) = mean(allFRs(Trials));
               FSEs(i,j)= std(allFRs(Trials)) / sqrt(sum(Trials));
           end
       end
    end

    rv = [FRs, FSEs, Values];
end
