function [ score ] = squeezeScore( d )
%squeezeScore gets a square matrix, which is a distance matrix
%   and calculates sum of the product of values and distances from the
%   diagonal

s = 0;
rawsum = 0;
for i = 1:size(d,1)
    for j = 1:size(d,1)
      s = s + d(i,j) * abs(i-j);  
      rawsum = rawsum + d(i, j);
    end
end

score = s / rawsum;
