function [e, t] = SpikesMulti(Vs, probe, threshold) 
 prePeak = 10;
 postPeak = 40;
 interSpikeGap = 30;
 ChSz = 10000;
 
 cc = 1;
 e = []; %zeros(length(p), prePeak + postPeak + 1);
 for chunk= 1 : 1 + length(Vs(probe,:))/ChSz
     if chunk*ChSz > length(Vs(probe,:))
         VsCh= Vs(:,1+(chunk-1)*ChSz:end);
     else
         VsCh= Vs(:,1+(chunk-1)*ChSz:chunk*ChSz);
     end
     disp(chunk)
    
     if (threshold < 0)
         p = find(VsCh(probe, :) < threshold);
         p = p(p>prePeak & p < length(VsCh(probe,:))-postPeak);
         
         for c = 1: length(p)
             if ((c == 1 ) || ((p(c) - p(c-1)) > interSpikeGap))
                 tp = find(VsCh(probe,p(c):p(c)+postPeak)==min(VsCh(probe,p(c)-prePeak:p(c)+postPeak)));
                     if p(c)+tp+postPeak < length(VsCh(probe,:)) 
                         e(cc, :, :) = VsCh(:,p(c)+tp-prePeak : p(c)+tp+postPeak);
                         t(cc) = p(c)+(chunk-1)*ChSz;
                         cc = cc + 1;
                         %if (rem(c,100) == 0)
                         %    disp([num2str(length(p)) , ' - ' , num2str(c) , ' - ' , num2str(cc)]);
                         %end
                     end
             end
         end
     else
         p = find(VsCh(probe, :) > threshold);
         p = p(p>prePeak & p < length(VsCh(probe,:))-postPeak);
         
         for c = 1: length(p)
             if ((c == 1 ) || ((p(c) - p(c-1)) > interSpikeGap))
                 tp = find(VsCh(probe,p(c):p(c)+postPeak)==max(VsCh(probe,p(c)-prePeak:p(c)+postPeak)));
                     if p(c)+tp+postPeak < length(VsCh(probe,:)) 
                         e(cc, :, :) = VsCh(:,p(c)+tp-prePeak : p(c)+tp+postPeak);
                         t(cc) = p(c)+(chunk-1)*ChSz;
                         cc = cc + 1;
                         %if (rem(c,100) == 0)
                         %    disp([num2str(length(p)) , ' - ' , num2str(c) , ' - ' , num2str(cc)]);
                         %end
                     end
             end
         end

         
     end
 end
 