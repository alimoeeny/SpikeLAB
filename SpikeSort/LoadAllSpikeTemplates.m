function [spikeTemplates] = LoadAllSpikeTemplates(spikeSize)


counter = 1;
spkTmplt = zeros(1,spikeSize);
spkTmplt(end/3) = 1;
spikeTemplates(counter,:) = spkTmplt;

counter = 2;
spkTmplt = zeros(1,spikeSize);
spkTmplt(end/2) = 1;
spikeTemplates(counter,:) = spkTmplt;

counter = 3;
spkTmplt = zeros(1,spikeSize);
spkTmplt(end*2/3) = 1;
spikeTemplates(counter,:) = spkTmplt;

counter = 4;
spkTmplt = zeros(1,spikeSize);
spkTmplt(end/3) = 1;
spkTmplt(end*2/3) = -1;
spikeTemplates(counter,:) = spkTmplt;

counter = 5;
spkTmplt = zeros(1,spikeSize);
spkTmplt(round(end/4)) = -1;
spkTmplt(round(end/2)) = 1;
spkTmplt(round(end*3/4)) = -1;
spikeTemplates(counter,:) = spkTmplt;

counter = 6;
spkTmplt = rand(1,spikeSize);
spikeTemplates(counter,:) = spkTmplt;


counter = 7;
spkTmplt = zeros(1,spikeSize);
spkTmplt = sin([1:spikeSize] ./ (spikeSize/ 6.28));
spikeTemplates(counter,:) = spkTmplt;

counter = 8;
spkTmplt = zeros(1,spikeSize);
spkTmplt = sin([1:spikeSize] ./ (spikeSize/ 3.141592));
spikeTemplates(counter,:) = spkTmplt;

counter = 9;
spkTmplt = zeros(1,spikeSize);
spkTmplt = sin([1:spikeSize] ./ (spikeSize/ 3.141592));
spkTmplt(6:end) = spkTmplt(1:end-5);
spkTmplt(1:5) = 0;
spkTmplt(end-5:end) = 0;
spikeTemplates(counter,:) = spkTmplt;

counter = 10;
spkTmplt = zeros(1,spikeSize);
spkTmplt = tan([1:spikeSize] ./ (spikeSize/ 6.282));
spkTmplt(6:end) = spkTmplt(1:end-5);
spkTmplt(1:5) = 0;
spkTmplt(end-5:end) = 0;
spikeTemplates(counter,:) = spkTmplt;


counter = 11;
spkTmplt = zeros(1,spikeSize);
spkTmplt = sin([1:spikeSize] ./ (spikeSize/ 9.42));
spikeTemplates(counter,:) = spkTmplt;

counter = 12;
spkTmplt = zeros(1,spikeSize);
spkTmplt = sin([1:spikeSize] ./ (spikeSize/ 15.7));
spkTmplt(1:5) = 0;
spkTmplt(end-5:end) = 0;
spikeTemplates(counter,:) = spkTmplt;

end

