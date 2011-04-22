function [e, t] = Spikes(V, threshold) 
 prePeak = 10;
 postPeak = 40;
 interSpikeGap = 30;
 
 if (threshold < 0)
    p = find(V < threshold);
    e = []; %zeros(length(p), prePeak + postPeak + 1);
    cc = 1;
    for c = 1: length(p)
        if ((c == 1 ) || ((p(c) - p(c-1)) > interSpikeGap))
            %disp([num2str(length(p)) , ' - ' , num2str(c) , ' - ' , num2str(cc)]);
            tp = find(V(p(c):p(c)+postPeak)==min(V(p(c)-prePeak:p(c)+postPeak)));
            e(cc, :) = V(p(c)+tp-prePeak : p(c)+tp+postPeak);
            t(cc) = p(c);
            cc = cc + 1;
            if (rem(c,100) == 0)
                disp([num2str(length(p)) , ' - ' , num2str(c) , ' - ' , num2str(cc)]);
            end
        end
    end
 else
     
 end
 