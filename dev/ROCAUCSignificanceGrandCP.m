function [rocauc1, pp, permutations] = ROCAUCSignificanceGrandCP(AA, BB)

if (length(AA) ~= length(BB))
    disp('invalid input for ROCAUCSignificanceGrandCP');
    return
end


aa = []; bb = [];
for i = 1:length(AA)
    aa = [aa; AA{i}];
    bb = [bb; BB{i}];
end
rocauc1 = ROCAUC(aa, bb);




numRep = 3000;
parfor rep = 1: numRep
    aap = []; bbp = [];
    for ip = 1:length(AA)
        aa = AA{ip};
        bb = BB{ip};
        t = [aa; bb];
        ri = randperm(length(t));
        aap = [aap; t(ri(1:length(aa)))];
        bbp = [bbp; t(ri(length(aa)+1:end))];
    end
       
    rocr(rep) = ROCAUC(aap, bbp);
end

rocr = sort(rocr);
[m, n] = find(rocr>=rocauc1);
if(size(n,2)>0)
    if(rocauc1> mean(rocr))
        pp = 100 - (n(1) * 100 / length(rocr));
    else
        pp = (n(1) * 100 / length(rocr));
    end
else
    pp = 0;
end

if nargout == 3
    permutations = rocr;
end



