function [ timestamps ] = detectEventsinRawVoltage( FullVV, vidx, max_spikeRate, direction)
% direction can be "up" (positive trigger) or "down" (negative voltage trigger)
% max_spikeRate is the maximum expected firing rate, determines the number of events returned 
FullV = cast(FullVV(vidx,:), 'double');
timestamps = [];
swatchLength = 30;
max_spikeCount = max_spikeRate * length(FullV) / 30000;%assuming 30K sample rate

dvdt = diff(FullV);
peakPoints = [];
peakPoints = find(dvdt==0);
dvdtXcross = dvdt(1:end-1) .* dvdt(2:end);
peakPoints = unique(sort([peakPoints find(dvdtXcross<0)]));
peakPoints = peakPoints(peakPoints>swatchLength);
peakPoints = peakPoints(peakPoints< (length(FullV) - swatchLength));
stdVV = std(FullV);
peakPoints = peakPoints(abs(FullV(peakPoints)) > stdVV);

peakvalues = FullV(peakPoints);

if(strcmpi(direction, 'up'))
    [~, peaksSortedIdx] = sort(peakvalues, 'descend');
else
    [~, peaksSortedIdx] = sort(peakvalues, 'ascend');
end

warning off
peakspointer = 1;
spikecounter = 0;
while ( (length(peakPoints)>peakspointer) && ...
        (length(timestamps)<max_spikeCount) && ...
        (spikecounter < length(peaksSortedIdx)))
    if((strcmpi(direction, 'up') && ...
        (FullV(peakPoints(peaksSortedIdx(peakspointer)))>= ...
            (max(FullV(peakPoints(peaksSortedIdx(peakspointer)) + [-swatchLength/2:swatchLength/2]))))) || ...
        (strcmpi(direction, 'down') && ...
        (FullV(peakPoints(peaksSortedIdx(peakspointer)))>= ...
            (min(FullV(peakPoints(peaksSortedIdx(peakspointer)) + [-swatchLength/2:swatchLength/2]))))))
        
            pp = peakPoints(peaksSortedIdx(peakspointer));
            
            sp = FullV(pp-swatchLength/2:pp+swatchLength/2);
            if(strcmpi(direction, 'up'))
                ppo = median(find(sp==max(sp), 1));
            else
                ppo = median(find(sp==min(sp), 1));
            end
        
            newtimestamp = pp - (swatchLength /2) + ppo;
        
            if ((sum(timestamps==newtimestamp)==0) && (sum(abs(timestamps - newtimestamp)<swatchLength/2)==0))
                timestamps(end+1) = newtimestamp;
                spikecounter = spikecounter + 1;
            end
    end
    peakspointer = peakspointer + 1;
end
warning on



disp(spikecounter);


end

