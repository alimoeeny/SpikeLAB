clear;
clc;

%AllPsychData = importdata('AllJBPsychDays.txt');
%AllPsychData = importdata('AllDaedalusPsychDays.txt');
%AllPsychData = {'/bgc/data/psych/dae/TwoCylF20Aug2010', '/bgc/data/psych/dae/TwoCylF23Aug2010'};
AllPsychData = {'/bgc/data/psych/dae/TwoCylF20Jul2011'};

%ExptType = 'bd';
ExptType = '';

rewardLevel = 0; % 0 for all, 2 for r1 3 for r2 4 for r3
if rewardLevel > 0 , disp(' O N L Y              S O M E           T R I A L S             I N C L U D E D !'); end

trialsLevel = 0; % 0 for all, 50 for first 50 of EACH BLOCK, 100 for first 100 of each block and ...
if trialsLevel > 0 , disp(' O N L Y              S O M E           T R I A L S             I N C L U D E D !'); end



%par
for iDay = 1: length(AllPsychData)
   disp (strcat(AllPsychData{iDay} ));
   if trialsLevel > 0
       ExptS = PsychMon(strcat(AllPsychData{iDay}),'getexpts');
       Expt = [];
       Expt.Stimvals = ExptS{1}.Stimvals;
       Expt.Trials = [];
       for es = 1: length(ExptS)
           if length(ExptS{es}.Trials) > trialsLevel
               for ti = 1: trialsLevel
                   Expt.Trials(end+1).dx    = ExptS{es}.Trials(ti).dx;
                   Expt.Trials(end).RespDir = ExptS{es}.Trials(ti).RespDir;
                   Expt.Trials(end).rw      = ExptS{es}.Trials(ti).rw;
               end
           end
       end
       
   else
       Expt = PsychMon(strcat(AllPsychData{iDay}),'getexpt'); 
   end
   if isempty(Expt)
        disp('Empty Experiment, no psych I guess!')
        continue
   end
   if rewardLevel > 0 ,
       rwValues = unique([Expt.Trials(:).rw]);
       if (length(rwValues)>=rewardLevel)
           Expt.Trials = Expt.Trials([Expt.Trials(:).rw]==rwValues(rewardLevel));
       else
           disp('N O T R I A L S in this R E W A R D L E V E L');
           continue;
       end
   end
   
   
   
   %Expt.Trials = Expt.Trials(1:100);
   disp (strcat(AllPsychData{iDay}, ' - Trials: ', num2str(length(Expt.Trials))));
   if strcmpi(ExptType, 'bd')
       dxValues = unique([Expt.Trials(:).bd]);
   else
       dxValues = unique([Expt.Trials(:).dx]);
   end
   dxValues = dxValues(dxValues<0.5 & dxValues>-0.5);
   tmpdxvals = [];
   for dxx = 1: length(dxValues)
        if strcmpi(ExptType, 'bd')
            if (sum([Expt.Trials(:).bd]==dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4 && sum([Expt.Trials(:).bd]== -dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4)
                tmpdxvals(end+1) = dxValues(dxx);
            end
        else
           if (sum([Expt.Trials(:).dx]==dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4 && sum([Expt.Trials(:).dx]== -dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4)
           tmpdxvals(end+1) = dxValues(dxx);
           end
        end
   end
   dxValues = tmpdxvals;     
   
   
   performances = zeros(length(dxValues),1);
   misses = zeros(length(dxValues),1);
   confintvs = zeros(length(dxValues), 1);
   for i = 1:length(dxValues),
      if strcmpi(ExptType, 'bd')
          performances(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
          misses(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == 0) / sum([Expt.Trials(:).bd] == dxValues(i));
          p = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == 1);
          n = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
      else
          %performances(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).score] == 1) / sum([Expt.Trials(:).dx] == dxValues(i) & ([Expt.Trials(:).score] == 1 | [Expt.Trials(:).score] == 0));
          performances(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == 1) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
          misses(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == 0) / sum([Expt.Trials(:).dx] == dxValues(i));
          p = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == 1);
          n = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
      end
      confintvs(i) = sqrt(p) * (1- p/ n);
   end
   
  
    
   psychdataset{iDay}.dxValues = dxValues;
   psychdataset{iDay}.performances = performances;
   psychdataset{iDay}.confintvs = confintvs;
   psychdataset{iDay}.or = [Expt.Stimvals.or];
   psychdataset{iDay}.missed = misses;
   try
       psychdataset{iDay}.slops = PsychSlop(Expt);
   catch 
       disp('something went wrong here');
   end
end

%% Graphics 

mycolors = 'rgbcmk'; mycolors = [mycolors mycolors mycolors mycolors mycolors mycolors mycolors mycolors mycolors mycolors mycolors];
figure,
for i = 1:length(AllPsychData) %length(AllPsychData)-5:length(AllPsychData)
    hold on,
    if(~isempty(psychdataset{i}))
        if(psychdataset{i}.or == 0) 
            mylinestyle = '-';
        else 
            mylinestyle = ':';
        end
        %plot(psychdataset{i}.dxValues, 100 - psychdataset{i}.performances, mycolors(i), 'LineWidth', 2, 'LineStyle', mylinestyle); 
        errorbar(psychdataset{i}.dxValues, 100 - psychdataset{i}.performances, psychdataset{i}.confintvs, mycolors(i), 'LineWidth', 2, 'LineStyle', mylinestyle); 
        %scatter(psychdataset{i}.dxValues, psychdataset{i}.performances, mycolors(i)); 
    end
end
refline(0);
refline(0, 50);
refline(0, 10);
%refline(0, 25);
refline(0, 90);
%refline(0, 75);
refline(0, 100);
%% 
alldxValues = [];
for i = 1:length(AllPsychData), 
    if(~isempty(psychdataset{i}))
        alldxValues = [alldxValues psychdataset{i}.dxValues]; 
    end
end
alldxValues = unique(alldxValues);

grandavg = zeros(length(alldxValues),2);
sarkhat =  zeros(length(alldxValues),2);
for i = 1:length(AllPsychData)
    for j = 1: length(alldxValues)
        if(~isempty(psychdataset{i}))
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
end

figure,
plot(grandavg./sarkhat);


%% Missed 
%for i = 1: 10, mean(psychdataset{i}.missed), end
figure,
bar(psychdataset{11}.missed);

