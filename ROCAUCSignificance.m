function [rocauc1, pp, permutations] = ROCAUCSignificance(a, b)

rocauc1 = ROCAUC(a, b);
numRep = 3000;
parfor rep = 1: numRep
    t = [a; b];
    ri = randperm(length(t));
    ap = t(ri(1:length(a)));
    bp = t(ri(length(a)+1:end));
    rocr(rep) = ROCAUC(ap, bp);
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