function [swatches, timestamps] = spikesOutResearch(trace, triggerTrace, thrshld, winSize)
swtchs = {};
tmstmps = {};
workerCount = 12;
presamples = 15;
postsamples = 15;

if(length(trace) ~= length(triggerTrace))
    disp('trace and triggertrace must be the same size');
    return
else
    %par
    for chunk = 0 : workerCount-1
        swtchsC = [];
        tmstmpsC = [];
        pointer = cast(winSize + winSize + (chunk * length(trace) / workerCount), 'int64');
        
        upperb = min(length(trace)- (2 * winSize), cast((chunk + 1) * length(trace) / workerCount, 'int64'));
        while pointer <  upperb
            if ( ((thrshld > 0) && (max(triggerTrace(pointer:pointer+winSize))> thrshld)) || ...
                 ((thrshld < 0) && (min(triggerTrace(pointer:pointer+winSize))< thrshld)))
                if(thrshld > 0)
                    mpoint = find(triggerTrace(pointer:pointer+winSize)==max(triggerTrace(pointer:pointer+winSize)), 1, 'first');
                    %mResearchPoint = find(trace(pointer-winSize:pointer+winSize)==max(trace(pointer-winSize:pointer+winSize)), 1, 'first');
                else
                    mpoint = find(triggerTrace(pointer:pointer+winSize)==min(triggerTrace(pointer:pointer+winSize)), 1, 'first');
                    %mResearchPoint = find(trace(pointer-winSize:pointer+winSize)==min(trace(pointer-winSize:pointer+winSize)), 1, 'first');
                end
                
                %mResearchPoint = mpoint + 15;
                %swtchsC(size(swtchsC,1)+1,:) = trace(pointer + mResearchPoint - presamples: pointer + mResearchPoint + postsamples);
                swtchsC(size(swtchsC,1)+1,:) = trace(pointer + mpoint - presamples: pointer + mpoint + postsamples);
                %tmstmpsC(size(tmstmpsC,1)+1) = mResearchPoint;
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