% Real Spike shape

clc
clear 


%% Load the data 

datapath =  '/sc/bgc5/Utah/lem/G002/';
%cd ~/Desktop/G001
goodchannels = [2 4 6 8 10 12 ] %[1:1:32]; %[1 3 5 7 9 11 13 15];
selectedchannel = 6;
i = 0;
for p = goodchannels%1:32
    i = i + 1;
    disp([datapath 'Expt1.p' num2str(p) 'FullV.mat']);
    load([datapath 'Expt1.p' num2str(p) 'FullV.mat'])
    v = FullV.V;
    full_voltage(i,:) = v(end-1000000:end); %(1:end/2);  %(v - mean(v)) ./ range(v);
    clear Full.V v
end

%%

%foldedSpace = reshape(full_voltage(:,1000001:2000000), size(full_voltage,1), 500, 2000);
foldedSpace = reshape(full_voltage(:,1:1000000), size(full_voltage,1), 500, 2000);
%foldedSpace = reshape(full_voltage(:,1:100000), size(full_voltage,1), 500, 200);
figure(29), plot(squeeze(mean(foldedSpace,3))')

%% calculate the average

averagedNoise = mean(full_voltage, 1);
figure(83), plot(averagedNoise);

%% calculate the gain and take the weighted average out
i = 0;
for p = goodchannels % 1:size(full_voltage,1)
    i = i + 1;
    g(i) = sum(full_voltage(i,:) .* averagedNoise) / sum(averagedNoise .* averagedNoise);
    disp([num2str(i) ' - gain is: ' num2str(g(i))]);
    filtered_voltage(i,:) = full_voltage(i,:) - (g(i) .* averagedNoise);
end

%% filter the unwanted frequencies out

freq_filtered_signal = eeg_filter(filtered_voltage(selectedchannel,:), 30000, 250, 10000, 2);
%hcfD = 250  / (30000 / 2);
%[BhpD, AhpD] = butter(4,hcfD, 'high');
%freq_filtered_signal = filtfilt(BhpD, AhpD, filtered_voltage(10,:));

figure(1819), clf, hold on, 
plot(full_voltage(selectedchannel,1:10000), 'r');
plot(averagedNoise(1:10000)', 'g', 'LineWidth', 4);
plot(freq_filtered_signal(1:10000)', 'k', 'LineWidth', 6);

%%

sliding_win_size = 150;
threshold = -140;
spikecounter = 1;
spikes = [];
slider = 100;
vcutoff = 2000;


for i = 1:length(filtered_voltage)-sliding_win_size
   if (i>slider)
    epoch = freq_filtered_signal(i:i+sliding_win_size); %filtered_voltage(10,i:i+sliding_win_size);%freq_filtered_signal(i:i+sliding_win_size);
    %epoch = (epoch - mean(epoch)) ./ range(epoch); %eeg_filter(epoch, 30000, 250, 10000, 2);
    %freq_filtered_epoch = freq_filtered_signal(i:i+sliding_win_size);%epoch; % eeg_filter(epoch, 30000, 250, 10000, 2);
    %if (max(freq_filtered_epoch)>abs(min(freq_filtered_epoch)))
       %local_peak = find(freq_filtered_epoch==max(freq_filtered_epoch));
       local_peak = find(epoch==min(epoch));
       if (local_peak(1)>1 && local_peak(1)<sliding_win_size*0.5 && epoch(local_peak(1))<threshold)
           spikes(spikecounter,1,:) = freq_filtered_signal(i+local_peak(1)-100:i+local_peak(1)+100) - freq_filtered_signal(i+local_peak(1)-20);
           spikes(spikecounter,2,:) = filtered_voltage(selectedchannel, i+local_peak(1)-100:i+local_peak(1)+100) - filtered_voltage(selectedchannel, i+local_peak(1)-20);
           spikes(spikecounter,3,:) = full_voltage(selectedchannel,i+local_peak(1)-100:i+local_peak(1)+100) - full_voltage(selectedchannel,i+local_peak(1)-20);
           
           spikeEnergies(spikecounter) = sum(diff(spikes(spikecounter,1,:)).^2);
           
           spikecounter = spikecounter + 1;
           slider = i + sliding_win_size;
           disp([num2str(i) ' - ' num2str(spikecounter)]);
       end
    %else
       
    %end
   end
end


%%
figure(999), clf, hold on,
hist(spikeEnergies,100);

%%
figure(9876), clf, hold on,
plot(squeeze(spikes(:,1,:))', 'b');
plot(squeeze(spikes(:,2,:))', 'g');
plot(squeeze(spikes(:,3,:))', 'r');

%%
figure(976), clf, hold on,
plot(squeeze(spikes(spikeEnergies>vcutoff,1,:))', 'b');
plot(squeeze(spikes(spikeEnergies<vcutoff,1,:))', 'c');
plot(squeeze(spikes(spikeEnergies>vcutoff,2,:))', 'g');
plot(squeeze(spikes(spikeEnergies<vcutoff,2,:))', 'y');
plot(squeeze(spikes(spikeEnergies>vcutoff,3,:))', 'r');
plot(squeeze(spikes(spikeEnergies<vcutoff,3,:))', 'm');


%%
figure(8765), clf, hold on,
plot(mean(squeeze(spikes(:,1,:)))', 'b');
plot(mean(squeeze(spikes(:,2,:)))', 'g');
plot(mean(squeeze(spikes(:,3,:)))', 'r');

%%
figure(875), clf, hold on,
plot(mean(squeeze(spikes(spikeEnergies>vcutoff,1,:)))', 'b');
plot(mean(squeeze(spikes(spikeEnergies<vcutoff,1,:)))', 'c');
plot(mean(squeeze(spikes(spikeEnergies>vcutoff,2,:)))', 'g');
plot(mean(squeeze(spikes(spikeEnergies<vcutoff,2,:)))', 'y');
plot(mean(squeeze(spikes(spikeEnergies>vcutoff,3,:)))', 'r');
plot(mean(squeeze(spikes(spikeEnergies<vcutoff,3,:)))', 'm');



%%

energycutoff=300000;
border = 240;
myspks = [];
spikeEnergies = [];

%load('/sc/bgc5/Utah/lem/G001/lemG001002.mat'); % -600
load('/sc/bgc5/Utah/lem/G002/lemG002001.mat');

peakVoltages = NEV.Data.Spikes.Waveform(:,13);
timeStamps = NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==goodchannels(selectedchannel));

timeStamps = timeStamps(timeStamps>2946739);

for i = 1:length(timeStamps)
    disp(i);
    spikeEnergies(i) = sum(diff(NEV.Data.Spikes.Waveform(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==goodchannels(selectedchannel)),:)).^2);
%     if  (peakVoltages(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==10))>-300)
%         spikeEnergies(i) = -spikeEnergies(i);
%     end 
    spks(i,:) = NEV.Data.Spikes.Waveform(find(NEV.Data.Spikes.TimeStamp==timeStamps(i) & NEV.Data.Spikes.Electrode==goodchannels(selectedchannel)),:); % - NEV.Data.Spikes.Waveform(i,10);
    myspks(i,1,:) = 1 * full_voltage(selectedchannel,timeStamps(i)-2946739-border:timeStamps(i)-2946739+border);
    myspks(i,2,:) = 1 * filtered_voltage(selectedchannel,timeStamps(i)-2946739-border:timeStamps(i)-2946739+border);
    myspks(i,3,:) = 1 * freq_filtered_signal(timeStamps(i)-border-2946739:timeStamps(i)-2946739+border);
    myspks(i,4,:) = 1 * averagedNoise(timeStamps(i)-border-2946739:timeStamps(i)+border-2946739);
end

figure(1983), clf, hold on, hist(spikeEnergies,100)

figure(1123), clf, hold on, 
plot(spks(spikeEnergies>energycutoff,:)', 'r', 'LineWidth', 2)
plot(spks(spikeEnergies<energycutoff,:)', 'k', 'LineWidth', 1)

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






%%

figure(6611), clf, hold on, 
plot(squeeze(myspks(spikeEnergies>energycutoff,1,:))', 'r', 'LineWidth', 2)
plot(squeeze(myspks(spikeEnergies<energycutoff,1,:))', 'k', 'LineWidth', 1)
reflinexy(border,200);

figure(6622), clf, hold on, 
plot(squeeze(myspks(spikeEnergies>energycutoff,2,:))', 'r', 'LineWidth', 2)
plot(squeeze(myspks(spikeEnergies<energycutoff,2,:))', 'k', 'LineWidth', 1)
reflinexy(border,200);


figure(6633), clf, hold on, 
plot(squeeze(myspks(spikeEnergies>energycutoff,3,:))', 'r', 'LineWidth', 2)
plot(squeeze(myspks(spikeEnergies<energycutoff,3,:))', 'k', 'LineWidth', 1)
reflinexy(border,200);


figure(6644), clf, hold on, 
plot(squeeze(myspks(spikeEnergies>energycutoff,4,:))', 'r', 'LineWidth', 2)
plot(squeeze(myspks(spikeEnergies<energycutoff,4,:))', 'k', 'LineWidth', 1)
reflinexy(border,200);



%%
% 
% close all, clear all
% 
% [zlp, plp, klp] = butter(3, 7500 / 15000  ,'low', 's');
% [soslp, glp] = zp2sos(zlp, plp, klp);
% Hd_LP = dfilt.df2tsos(soslp, glp);
% 
% [zhp, php, khp] = butter(1, 30 / 15000.0 ,'high', 's');
% [soshp, ghp] = zp2sos(zhp, php, khp);
% Hd_HP = dfilt.df2tsos(soshp, ghp);
% 
% [zDhp, pDhp, kDhp] = butter(4, 250.0 / 15000.0 ,'high');
% [sosDhp, gDhp] = zp2sos(zDhp, pDhp, kDhp);
% Hd_D_HP = dfilt.df2tsos(sosDhp, gDhp);
% 
% [z, p, k] = butter(3, [30 7500.0] / 15000.0 ,'bandpass');
% [sos, g] = zp2sos(z, p, k);
% Hd = dfilt.df2tsos(sos, g);
% %h = fvtool(Hd, 'FrequencyScale', 'log');
% %set(h, 'Analysis', 'freq')
% 
% %h = fvtool(Hd, Hd_LP, Hd_HP, Hd_D_HP, 'FrequencyScale', 'log');
% %set(h, 'Analysis', 'freq')
% 
% h = fvtool(Hd_LP, 'FrequencyScale', 'log');
% 


%% blackrock 

lcf = 0.3  / (30000 / 2);
hcf = 7500 / (30000 / 2);

% [B, A] = butter(3,[lcf, hcf]);
% blackrockBand = filtfilt(B, A, filtered_voltage(10,:));

[Bhp, Ahp] = butter(1,lcf, 'high');
[Blp, Alp] = butter(3,hcf, 'low');
hcfD = 250  / (30000 / 2);
[BhpD, AhpD] = butter(4,hcfD, 'high');
blackrock = filtfilt(BhpD, AhpD, filtfilt(Bhp, Ahp, filtfilt(Blp, Alp, filtered_voltage(10,:))));



%% 
winsize = 48;
padding = 20;
spikecounter = 1;

spikes = [];
rawspikes = [];
spikesBlkRck = [];
%spikesBlkRckBND = [];

pk_i = 0;
figure(1323), clf, hold on
thrshld = 58.4; %95; %250
energyThreshold = 2350;

signal = freq_fltered_signal; %blackrock; %freq_fltered_signal;

peaks = find(signal>thrshld);
slopes = diff(signal);

for i = 1:length(peaks)
    disp(i);
    
    epoch = signal(peaks(i)-10-padding:peaks(i)+winsize+padding);
    
    if (mod(spikecounter, 100)==0)
        debug=1;
    end
    
    if (peaks(i)>pk_i+30.)
        pk = peaks(i);
        pk_i = peaks(i);
        spikes(spikecounter,:) = epoch;
        spikeEnergies(spikecounter) = sum(diff(spikes(spikecounter,:)).^2);
        rawspikes(spikecounter,:) = full_voltage(10, peaks(i)-10-padding:peaks(i)+winsize+padding);
        spikesBlkRck(spikecounter,:) = blackrock(peaks(i)-10-padding:peaks(i)+winsize+padding);
               
        if(spikeEnergies(spikecounter)>energyThreshold)
            plot(spikes(spikecounter,:), 'b');
            plot(rawspikes(spikecounter,:), 'r');
            plot(spikesBlkRck(spikecounter,:), 'g');
        end
        spikecounter = spikecounter + 1;
    end
end

%%
figure(8346), hist(spikeEnergies, 100);

%%
figure(1818), clf, hold on
%plot(mean(spikes,1)), 
plot(mean(spikes(spikeEnergies>energyThreshold,:),1), 'm', 'LineWidth', 2),
plot(mean(spikes(spikeEnergies<energyThreshold,:),1), 'k', 'LineWidth', 2),
plot(mean(rawspikes(spikeEnergies>energyThreshold,:),1), 'm--', 'LineWidth', 2),
plot(mean(rawspikes(spikeEnergies<energyThreshold,:),1), 'k--', 'LineWidth', 2),

plot(mean(spikesBlkRck(spikeEnergies>energyThreshold,:),1), 'r', 'LineWidth', 2),
plot(mean(spikesBlkRck(spikeEnergies<energyThreshold,:),1), 'g', 'LineWidth', 2),

reflinexy(padding+10,1000);
reflinexy(padding+40,1000);
reflinexy(padding+70,1000);
reflinexy(padding+100,1000);


%% 

spikespace = [];

[COEFF,SCORE] = princomp(spikes);
for s = 1:size(spikes, 1)
   spikespace(s,1) =  var(spikes(s,:)) ./ sum(abs(spikes(s,:)));
   spikespace(s,2) = var(spikes(s,:));
end
figure, scatter(spikespace(:,1), spikespace(:,2))



%%

%Hhp = fdesign.highpass('Fst,Fp,Ast,Ap',10,15,60,1,30000);



%%
figure, 
r = [1:length(full_voltage)]; plot (full_voltage(10, r)), hold on, plot(filtered_voltage(10, r), 'r'), hold on, plot(averagedNoise(r), 'g');
r = [1:1000];  plot (full_voltage(10, r)), hold on, plot(filtered_voltage(10, r), 'r'), hold on, plot(averagedNoise(r), 'g');
%plot (full_voltage(10, 1:10000)), hold on, plot(filtered_voltage(10, 1:10000), 'r');


%%
figure, periodogram(averagedNoise(1:end),[],2^16,30E3);
hold on, periodogram(filtered_voltage(10,1:end),[],2^16,30E3);
hold on, periodogram(full_voltage(10,1:end),[],2^16,30E3);

%%

window_size = 40; % slide a window of 48 samples over the data and detect peaks

hp_kernel = ones(1800,1);
%hp_kernel = gausswin(600,2.5)-gausswin(600,1.25);
%hp_kernel = hp_kernel(301:end);

filtered_voltage = conv(full_voltage, hp_kernel);

plot (full_voltage(1:10000)), hold on, plot(filtered_voltage(1:10000), 'r')

figure, spectrogram(full_voltage(1:1000), 30,15, 300, 30E3)



