clear,clc;

% AllDaedalusPsychDays = importdata('AllDaedalusPsychDays.txt');
load('AllDaeBDIDPsychDays.mat');
AllDaedalusPsychDays = AllDaeBDIDPsychDays;
clear AllDaeBDIDPsychDays;

%psychdataset = cell(1,length(AllDaedalusPsychDays));

psychPath = '/sc/bgc/data/psych/dae/';

Daycntr = 0;
%par
for iDay = 1: length(AllDaedalusPsychDays)
    disp(AllDaedalusPsychDays(iDay));

    ExptS = PsychMon([psychPath, AllDaedalusPsychDays{iDay}],'getexpts');
   
    for e = 1: length(ExptS)
        Expt = ExptS{e};
        if (length(Expt.Trials)>70)
            if ~((isfield(Expt.Trials, 'bd')) && (isfield(Expt.Trials, 'Id')))
               disp( 'NO bd and Id HERE!');
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

                 or = Expt.Stimvals.or; % ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
                 if (or==90)
                     c1 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                     c2 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                     c3 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                     c4 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                 elseif (or==-90)
                     c1 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                     c2 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                     c3 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                     c4 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);             
                 elseif (or==0)
                     c1 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                     c2 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);
                     c3 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                     c4 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);             
                 elseif (or==180)
                     c1 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] > 0);
                     c2 = ([Expt.Trials(:).Id] < 0) & ([Expt.Trials(:).RespDir] < 0);
                     c3 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] > 0);
                     c4 = ([Expt.Trials(:).Id] > 0) & ([Expt.Trials(:).RespDir] < 0);             
                 else 
                     c1 = 0; c2 = 0; c3 = 0; c4 = 0;
                     disp(['ORIENTATION IS :  ', num2str(or), '  WWHHAATT CCAANN II DDOO == == == == == == == =='])
                 end
                 bdCrossTalk5 = 0.5 * ( ...
                                    ((sum(c1) - sum(c2)) / (sum(c1) + sum(c2))) + ...
                                    ((sum(c4) - sum(c3)) / (sum(c4) + sum(c3))) );


                                
                slp = PsychSlop(Expt, 'bd', 'BDID');
                                
                Daycntr = Daycntr + 1;
                psychdataset{Daycntr}.dxValues = dxValues;
                psychdataset{Daycntr}.performances = performances;
                psychdataset{Daycntr}.or = [Expt.Stimvals.or];
                psychdataset{Daycntr}.missed = misses;   
                psychdataset{Daycntr}.slope = slp.fit(2);
                psychdataset{Daycntr}.bdCrossTalk5 = bdCrossTalk5;
            end
        end
    end
end


%% Graphics
mycolors = 'rgbcmk'; mycolors = [mycolors mycolors mycolors mycolors mycolors mycolors mycolors];

for i = 1:length(psychdataset)
    if ~isempty(psychdataset{i})
        or(i) = psychdataset{i}.or;
        if ((psychdataset{i}.slope > 100) | (psychdataset{i}.slope < -100))
            slope(i) = -0.987654321;
        else
            slope(i) = psychdataset{i}.slope;
        end
        bdCrossTalk5(i) = psychdataset{i}.bdCrossTalk5;
    else 
        or(i) = -500;
        slope(i) = -500;
        bdCrossTalk5 = -500;
    end
end

%% 
figure,
for i = 1:length(psychdataset)
    hold on,
    if ~isempty(psychdataset{i})
        if(psychdataset{i}.or == 0) 
            mylinestyle = '-';
        else 
            mylinestyle = ':';
        end
        plot(psychdataset{i}.dxValues, psychdataset{i}.performances, mycolors(i), 'LineWidth', 2, 'LineStyle', mylinestyle); 
        %scatter(psychdataset{i}.dxValues, psychdataset{i}.performances, mycolors(i)); 
    end
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

