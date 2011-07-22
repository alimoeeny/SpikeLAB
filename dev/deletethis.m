
for iN = 1: length(AllNeurons)
    [a,x, y, z] = ttest(PSTHs{1,iN}{3}(1,:) ,PSTHs{1,iN}{3}(2,:), 0.001);
    b = mean(PSTHs{1,iN}{3}(1,:) - PSTHs{1,iN}{3}(2,:)) / mean(PSTHs{1,iN}{3}(2,:));
    p = anova1([PSTHs{1,iN}{3}(1,:); PSTHs{1,iN}{3}(2,:); PSTHs{1,iN}{3}(3,:)]', [], 'off');
    
    conditions = PSTHs{1,iN}{2};
    
    d =  mean([mean(PSTHs{1,iN}{1}(conditions(3,:),:),2) ; mean(PSTHs{1,iN}{1}(conditions(4,:),:),2)]) - ...
         mean([mean(PSTHs{1,iN}{1}(conditions(5,:),:),2) ; mean(PSTHs{1,iN}{1}(conditions(6,:),:),2)]) ;
    e =  ( mean(mean(PSTHs{1,iN}{1}(conditions(1,:),:),2)) - mean(mean(PSTHs{1,iN}{1}(conditions(2,:),:),2)) );
    
    %d = mean( (PSTHs{1,iN}{3}(3,:) + PSTHs{1,iN}{3}(4,:)) - (PSTHs{1,iN}{3}(5,:) + PSTHs{1,iN}{3}(6,:)) ) * 0.5 ;
    %r = mean(  PSTHs{1,iN}{3}(1,:) - PSTHs{1,iN}{3}(2,:)) / d;
    
    disp([abs(TI(iN))>0.1,a, b, p, e, e / d]);
    aa = a;
    bb(iN) = b;
    ee(iN) = e;
    dd(iN) = d;
end