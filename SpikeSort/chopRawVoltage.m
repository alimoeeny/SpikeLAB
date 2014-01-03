function [swatches] = chopRawVoltage(FullVV, vidx, timestamps)
 
swatches = [];
swatchLength = 30;

for it =1:length(timestamps)
    swatches(it, :) = FullVV(vidx, timestamps(it)+[-swatchLength/2:swatchLength/2]);
end

end

