function [ sixteen ] = eight2sixteen( eight )
for i = 2:2:length(eight)
    sixteen(i/2) =  double(eight(i-1)) * 256.0 + double(eight(i));
end

end

