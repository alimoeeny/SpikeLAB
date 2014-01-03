function [swatches, timestamps] = spikesOut(trace, triggerTrace, thrshld, winSize)
swtchs = {};
tmstmps = {};
workerCount = 12;
presamples = 15;
postsamples = 15;

if(size(trace) ~= size(triggerTrace))
    disp('trace and triggertrace must be the same size');
    return
else
    %par
    for chunk = 0 : workerCount-1
        swtchsC = [];
        tmstmpsC = [];
        pointer = cast(winSize + (chunk * length(trace) / workerCount), 'int64');
        
        upperb = min(length(trace)-winSize, cast((chunk + 1) * length(trace) / workerCount, 'int64'));
        while pointer <  upperb
            if ( ((thrshld > 0) && (max(triggerTrace(pointer:pointer+winSize))> thrshld)) || ...
                 ((thrshld < 0) && (min(triggerTrace(pointer:pointer+winSize))< thrshld)))
             if(thrshld > 0)   
                 mpoint = find(triggerTrace(pointer:pointer+winSize)==max(triggerTrace(pointer:pointer+winSize)), 1, 'first');
             else
                 mpoint = find(triggerTrace(pointer:pointer+winSize)==min(triggerTrace(pointer:pointer+winSize)), 1, 'first');
             end
                swtchsC(size(swtchsC,1)+1,:) = trace(pointer + mpoint - presamples: pointer + mpoint + postsamples);
                tmstmpsC(size(tmstmpsC,1)+1) = mpoint;
                pointer = pointer + winSize + 15;
                %disp(['chunk: ', num2str(chunk)]);
            else
                pointer = pointer + winSize;
            end
        end
        
        swtchs{chunk+1} = swtchsC;
        tmstmps{chunk+1}= tmstmpsC;
    end
    
    swatches = [];
    timestamps = [];
    for chunk = 1: workerCount
        if(~isempty(swtchs{chunk}))
            swatches = [swatches; swtchs{chunk}];
            timestamps=[timestamps; tmstmps{chunk}'];
        end
    end
end
end