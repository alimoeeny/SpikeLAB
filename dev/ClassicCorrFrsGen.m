clear, clc,

PoolCount = 5;
nn = 100;%1000 
Nn = PoolCount * nn;
Tn = 100;
PopmFr = 1;


disp('Let''s do this!');
PoolCountHalf = 4;
corrWithin = 0.2;
corrMinBetween = 0.1;
CorrMat = zeros(Nn, Nn) + corrMinBetween;
for pool = 1: PoolCount
    for offset = 1: PoolCount
        if (((pool > 4) || (offset > 4)) && (pool~= offset))
            CorrMat((pool-1)*nn +1 : pool*nn, (offset-1)*nn + 1 : (offset)*nn) = corrMinBetween;
        else
            CorrMat((pool-1)*nn +1 : pool*nn, (offset-1)*nn + 1 : (offset)*nn) = corrWithin - (1.0*abs(-pool + offset)/(PoolCountHalf-1) * (corrWithin - corrMinBetween));
        end
    end
end

for i = 1: Nn
        CorrMat(i, i) = 1;
end       

disp('CorrMat READY!'); 
Frs = normrnd(PopmFr, PopmFr, [Tn, Nn]);
disp('Frs READY!'); 
CorrMatSq = sqrtm(CorrMat);
disp('CorrMatSq READY!'); 
Frs = Frs * CorrMatSq;
disp('Saving ...'); 
save('/Users/ali/Desktop/FrsPyClassGen.mat', 'Frs');

%% Sanity!
figure(345), imagesc(CorrMat);
figure(357), imagesc(Frs);
figure(444), imagesc(CorrMatSq);
cimsc = parcorr(Frs');
figure(789), imagesc(cimsc);