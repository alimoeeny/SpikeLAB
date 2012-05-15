
energycutoff=100000;
border = 240;
myspks = [];
spikeEnergies = [];

%load('/sc/bgc5/Utah/lem/G001/lemG001002.mat'); % -600
load('/sc/bgc5/Utah/lem/G002/lemG002001.mat');

peakVoltages = NEV.Data.Spikes.Waveform(:,13);
%probe = 1;
for probe = 1:24
timeStamps = NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==probe);


for i = 1:length(timeStamps)
    disp(i);
    spikeEnergies(i) = sum(diff(NEV.Data.Spikes.Waveform(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==probe),:)).^2);
%     if  (peakVoltages(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==10))>-300)
%         spikeEnergies(i) = -spikeEnergies(i);
%     end 
    spks(i,:) = NEV.Data.Spikes.Waveform(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==probe),:); % - NEV.Data.Spikes.Waveform(i,10);
%     myspks(i,1,:) = 1 * full_voltage(probe,timeStamps(i)-border:timeStamps(i)+border);
%     myspks(i,2,:) = 1 * filtered_voltage(probe,timeStamps(i)-border:timeStamps(i)+border);
%     myspks(i,3,:) = 1 * freq_filtered_signal(timeStamps(i)-border:timeStamps(i)+border);
%     myspks(i,4,:) = 1 * averagedNoise(timeStamps(i)-border:timeStamps(i)+border);
end

figure(1982+probe), clf, hold on, hist(spikeEnergies,100)

figure(1122+probe), clf, hold on, 
plot(spks(spikeEnergies>energycutoff,:)', 'r', 'LineWidth', 2)
plot(spks(spikeEnergies<energycutoff,:)', 'k', 'LineWidth', 1)

end

figure(3311), clf, hold on, 
plot(squeeze(mean(myspks(spikeEnergies>energycutoff,1,:),1))', 'r', 'LineWidth', 2)
plot(squeeze(mean(myspks(spikeEnergies<energycutoff,1,:),1))', 'k', 'LineWidth', 1)
reflinexy(border,200);

figure(3322), clf, hold on, 
plot(squeeze(mean(myspks(spikeEnergies>energycutoff,2,:),1))', 'r', 'LineWidth', 2)
plot(squeeze(mean(myspks(spikeEnergies<energycutoff,2,:),1))', 'k', 'LineWidth', 1)
reflinexy(border,200);

figure(3333), clf, hold on, 
plot(squeeze(mean(myspks(spikeEnergies>energycutoff,3,:),1))', 'r', 'LineWidth', 2)
plot(squeeze(mean(myspks(spikeEnergies<energycutoff,3,:),1))', 'k', 'LineWidth', 1)
reflinexy(border,200);

figure(3344), clf, hold on, 
plot(squeeze(mean(myspks(spikeEnergies>energycutoff,4,:),1))', 'r', 'LineWidth', 2)
plot(squeeze(mean(myspks(spikeEnergies<energycutoff,4,:),1))', 'k', 'LineWidth', 1)
reflinexy(border,200);
