function ss = squeezeOneStep(spc, step)
figure(181+step), clf, subplot(1,2,1); 
imagesc(tril(distancesMatrix(spc)));
s = spc([step:end+1-step], :); 
d = distancesMatrix(s);
[m, midx] = max(d(:));
[a, b] = ind2sub(size(d), midx);
if(a>b), tmp = a; a = b; b = tmp; end
if (a==1)
    tmpspace = [s([1:b-1 b+1:end],:); s(b,:)];
else
    tmpspace = [s(a,:); s([1:a-1 a+1:b-1 b+1:end],:); s(b,:)];
end

spc([step:end+1-step], :) = tmpspace;
subplot(1,2,2), imagesc(tril(distancesMatrix(spc)));
ss = spc;
end
