function d = distancesMatrix(s)
d = zeros(size(s,1), size(s,1));
for i = 1: size(s,1)
    for j = 1:size(s,1)
        d(i,j) = sqrt((s(i,1) - s(j,1))^2 + (s(i,2) - s(j,2))^2);
    end
end
end