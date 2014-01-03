function r = dismat(uselessinput)
space = [1 3; 2 1; 3 1; 12 8; 12 11; 10 9; 1.5 2; 1 2; 1.75 1.5; 1.65 2.1; 10.5 11; 2 3];
tmpspace = [];
newdistances = [];

figure(1122), clf, hold on,
s1 = squeezeOneStep(space,1);
s2 = squeezeOneStep(s1, 2);
s3 = squeezeOneStep(s2, 3);
s4 = squeezeOneStep(s3, 4);
s5 = squeezeOneStep(s4, 5);

end





