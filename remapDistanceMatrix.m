function [remapedDM] = remapDistanceMatrix(dm, map)

remapedDM = zeros(size(dm));
for i = 1:size(dm,1)
    for j  = 1:size(dm,2)
        remapedDM(i,j) = dm(map(i), map(j));
    end
end

end