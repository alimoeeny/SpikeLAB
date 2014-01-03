function [scores, permutations ] = forceSqueezeDistanceMatrix( originalDistanceMatrix )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%space = [1 3; 2 1; 3 1; 12 8; 12 11; 10 9; 1.5 2; 1 2; 1.75 1.5; 1.65 2.1; 10.5 11; 2 3];
%space = [1 3; 2 1; 5.5 4.5; 12 8; 5 6; 10 9; 1.5 2; 1 2; 1.75 1.5];
%space = [1 3;  3 1; 12 11;  1.5 2; 1.75 1.5; 10.5 11];
%figure(1112), scatter(space(:,1), space(:,2), 'r', 'filled')
%p = perms(1:size(space,1));
%initialDistances = distancesMatrix(space);
%initialScore = squeezeScore(initialDistances);

p = perms(1:size(originalDistanceMatrix,1));
scores = [];

for ip = 1:size(p,1)
    remapedDM = remapDistanceMatrix(originalDistanceMatrix, p(ip,:));
    scores(ip) = squeezeScore(remapedDM);
end

permutations = p;
% maxidx = find(scores==max(scores));
% figure(959), imagesc(remapDistanceMatrix(initialDistances, p(maxidx,:)));
% minidx = find(scores==min(scores));
% figure(969), imagesc(remapDistanceMatrix(initialDistances, p(minidx,:)));
end
