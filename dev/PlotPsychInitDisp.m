clear,clc;

% AllDaedalusPsychDays = importdata('AllDaedalusPsychDays.txt');
load('AllDaeBDIDPsychDays.mat');
AllDaedalusPsychDays = AllDaeBDIDPsychDays;
clear AllDaeBDIDPsychDays;

%psychdataset = cell(1,length(AllDaedalusPsychDays));

psychPath = '/sc/bgc/data/psych/dae/';

%par
parfor iDay = 1: length(AllDaedalusPsychDays)
   
    disp(AllDaedalusPsychDays(iDay));

    Expt = PsychMon([psychPath, AllDaedalusPsychDays{iDay}],'getexpt');
   
    if ~(isfield(Expt.Trials, 'bd'))
       disp( 'NO bd HERE!');
        continue;
    else
        dxValues = unique([Expt.Trials(:).bd]);
    
        performances = zeros(length(dxValues),1);
        misses = zeros(length(dxValues),1);
        for i = 1:length(dxValues),
        %performances(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).score] == 1) / sum([Expt.Trials(:).dx] == dxValues(i) & ([Expt.Trials(:).score] == 1 | [Expt.Trials(:).score] == 0));
        %performances(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
        performances(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
        %misses(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == 0) / sum([Expt.Trials(:).dx] == dxValues(i));
        misses(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == 0) / sum([Expt.Trials(:).bd] == dxValues(i));
        end
   
  
   
        psychdataset{iDay}.dxValues = dxValues;
        psychdataset{iDay}.performances = performances;
        psychdataset{iDay}.or = [Expt.Stimvals.or];
        psychdataset{iDay}.missed = misses;   
    end
end

%% Graphics
mycolors = 'rgbcmk'; mycolors = [mycolors mycolors mycolors mycolors mycolors mycolors mycolors];
figure,
for i = 1:length(psychdataset)
    hold on,
    if(psychdataset{i}.or == 0) 
        mylinestyle = '-';
    else 
        mylinestyle = ':';
    end
    plot(psychdataset{i}.dxValues, psychdataset{i}.performances, mycolors(i), 'LineWidth', 2, 'LineStyle', mylinestyle); 
    %scatter(psychdataset{i}.dxValues, psychdataset{i}.performances, mycolors(i)); 
end


alldxValues = [];
for i = 1:length(psychdataset), alldxValues = [alldxValues psychdataset{i}.dxValues]; end
alldxValues = unique(alldxValues);

grandavg = zeros(length(alldxValues),2);
sarkhat =  zeros(length(alldxValues),2);
for i = 1:length(AllDaedalusPsychDays)
    for j = 1: length(alldxValues)
        tempvalue = psychdataset{i}.performances(psychdataset{i}.dxValues == alldxValues(j));
        if(~isnan(tempvalue)) 
            if(psychdataset{i}.or == 0)
                grandavg(j,1) = grandavg(j,1) + psychdataset{i}.performances(psychdataset{i}.dxValues == alldxValues(j)); 
                sarkhat(j,1) = sarkhat(j,1) + 1;
            else
                grandavg(j,2) = grandavg(j,2) + psychdataset{i}.performances(psychdataset{i}.dxValues == alldxValues(j)); 
                sarkhat(j,2) = sarkhat(j,2) + 1;
            end
        end
    end
end

figure,
plot(grandavg./sarkhat);


%% Missed 
for i = 1: 10, mean(psychdataset{i}.missed), end
figure,
bar(psychdataset{11}.missed);

