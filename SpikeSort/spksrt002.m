
template1 = zeros(1,30); template1([7 8]) = -1; template1([15 16]) = 1;
template2 = zeros(1,30); template2([7 8]) = 1; template2([15 16]) = -1;
template3 = zeros(1,30); template3([5 6]) = 1; template3([9 10]) = -1; template3([15 16]) = 1;
template4 = zeros(1,30); template4([1 2]) = 1; template4([5 6]) = -1; template4([9 10]) = 1; template4([13 14]) = -1;

%load /Volumes/bgc5/smr/lem/M221/Expt10FullV.mat
%load /data/Expt6FullV.mat
V12 = cast(FullV.V(12,:), 'double') ;
V12n = V12 ./ (2 * std(V12));

V12nT1 = conv(V12n, template1);
V12nT1 = V12nT1 ./ (2 * std(V12nT1));

V12nT2 = conv(V12n, template2);
V12nT2 = V12nT2 ./ (2 * std(V12nT2));

V12nT3 = conv(V12n, template3);
V12nT3 = V12nT3 ./ (2 * std(V12nT3));

V12nT4 = conv(V12n, template4);
V12nT4 = V12nT4 ./ (2 * std(V12nT4));


[s, t] = spikesOut(V12nT4, V12nT4, 3, 30);

figure, plot(mean(s))
hold on, plot(s')