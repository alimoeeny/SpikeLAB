clc;
clear;
DataPath = GetDataPath('server');


ShowIndividualPlots = 0; % 0 or 1
SaveIndividualPlots = 0; % 0 or 1
[AllNeurons, FileType, StimulusType] = loadAllNeurons4('TWO');

% Psych
% AllPsychData = importdata('AllDaedalusPsychDays.txt');
% AllNeurons = AllPsychData;
% clear AllPsychdata;
% FileType = 'psych';
% StimulusType = 'cylinder';


% ABD
%  load('/Users/ali/DropBox/Projects/BCode/AllABDNeurons.mat');
%  AllNeurons = AllABDNeurons;
%  clear AllABDNeurons;
%  FileType = 'ABD';
%  StimulusType = 'cylinder';

% % DID
% load('../AllDIDNeurons.mat');
% AllNeurons = AllDIDNeurons;
% clear AllDIDNeurons;
% FileType = 'DID';
% StimulusType = 'cylinder';

% load('/Users/ali/DropBox/Projects/BCode/AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWONeurons;
% FileType = 'TWO';
% StimulusType = 'cylinder';

% TWO
% load('/Users/ali/DropBox/Projects/BCode/AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWOPsychDays;
% FileType = 'TWO';
% StimulusType = 'cylinder';

% PSYCH TWO
% load('/Users/ali/DropBox/Projects/BCode/AllTWOPsychDays.mat');
% AllNeurons = AllTWOPsychDays;
% clear AllTWOPsychDays;
% FileType = 'TWOPsych';
% StimulusType = 'cylinder';

% BD DID psych
% AllNeurons = {'/bgc/data/psych/dae/TwoCylF07Jun2010'};
% FileType = 'BDpsych';
% StimulusType = 'cylinder';

% load('/Users/ali/DropBox/Projects/BCode/AllRIDNeurons.mat');
% AllNeurons = AllRIDNeurons;
% clear AllRIDNeurons;
% FileType = 'RID';
% StimulusType = 'rds';


StartTime = 5500; 
FinishTime = 0; 

%SelectedNeurons = 1:length(AllNeurons); 
%par
for iN= 1:length(AllNeurons), %iN= [20: 25] %length(AllNeurons)] %[1:20] % 1:length(AllNeurons), 
    NeuronName = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronName); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));

%    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat'); 
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);

    TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, [], DataPath);
    if(abs(TI(iN))>0.05) 
    
        
        if(~isempty(strfind(FileType, 'psych')) || strcmp(FileType, 'TWOPsych'))
            Expt = PsychMon(strcat(AllNeurons{iN}),'getexpt');
        else
            Neuron = load(filepath);
            Expt = Neuron.Expt;
        end
        if isempty(Expt)
            disp('Empty Experiment, no psych I guess!')
            continue
        end
        if(~isfield(Expt.Trials, 'Id') && (strcmp(FileType, 'ABD') || strcmp(FileType, 'DID') ))
            disp('Not an idisp expt');
            continue;
        end
        if strcmpi(FileType, 'BDpsych')
            dxValues = unique([Expt.Trials(:).bd]);
        else
            dxValues = unique([Expt.Trials(:).dx]);
        end
        if (strcmp(FileType, 'psych'))
            dxValues = dxValues(dxValues<0.5 & dxValues>-0.5);
            tmpdxvals = [];
            for dxx = 1: length(dxValues)
                if (sum([Expt.Trials(:).dx]==dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4 && sum([Expt.Trials(:).dx]== -dxValues(dxx) & [Expt.Trials(:).RespDir] ~= 0) > 4)
                    tmpdxvals(end+1) = dxValues(dxx);
                end
            end
            dxValues = tmpdxvals;
            if(length(dxValues)<3 || sum(dxValues==0)<1)
                disp('Too few dx values or no dx==0');
                continue;
            end
        end
        if strcmpi(FileType, 'BDpsych')
            if(mean([Expt.Trials([Expt.Trials(:).bd]>0).RespDir])>0)
                ResponseToPositive = 1;
                ResponseToNegative = -1;
            else
                ResponseToPositive = -1;
                ResponseToNegative = 1;
            end
        else
            if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
                ResponseToPositive = 1;
                ResponseToNegative = -1;
            else
                ResponseToPositive = -1;
                ResponseToNegative = 1;
            end
        end
        dxs{iN} = dxValues;
        if(iseven(size(dxValues,2))), disp(strcat('It is even, it is odd!', ' - ' , num2str(iN))), end
        if (strcmp(FileType, 'TWO') || strcmpi(FileType, 'BDpsych'))
            ids{iN} = unique([Expt.Trials(:).bd]);
        else
            ids{iN} = unique([Expt.Trials(:).Id]);
        end
        Grand = [];
        dxss = [];
        trials = [];
        
        if strcmp(FileType, 'DID')
            disp([NeuronNumber, max(dxValues) , max([Expt.Trials(:).Id])]);
        end
        
        i=1;
        if strcmpi(FileType, 'BDpsych')
            Grand(1) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            dxss(1) = dxValues(i);
            trials(1) = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
        else if strcmp(FileType, 'TWO')
                Grand(1) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                Grand(11)= 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == -dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == -dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                if(NeuronNumber == 324 || NeuronNumber == 330)
                    Grand(11)= 100 * sum([Expt.Trials(:).dx] == dxValues(i+1) & [Expt.Trials(:).bd] == -dxValues(i+1) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i+1) & [Expt.Trials(:).bd] == -dxValues(i+1) & [Expt.Trials(:).RespDir] ~= 0);
                end
                dxss(1) = dxValues(i);
                trials(1) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            else
                Grand(1) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                dxss(1) = dxValues(i);
                trials(1) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            end
        end
        g= []; d= []; t=[];
        for i = [2:1:size(dxValues,2)/2]
            if strcmpi(FileType, 'BDpsych')
                g(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                d(i) = dxValues(i);
                t(i) = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            else if strcmp(FileType, 'TWO')
                    g(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                    d(i) = dxValues(i);
                    t(i) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                else
                    g(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                    d(i) = dxValues(i);
                    t(i) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                end
            end
        end
        Grand(2) = mean(g([2:1:size(dxValues,2)/2]));
        dxss(2) = mean(d([2:1:size(dxValues,2)/2]));
        trials(2) = sum(t([2:1:size(dxValues,2)/2]));
        
        
        if strcmpi(FileType, 'BDpsych')
            Grand(3) = 100 * sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0);
            Grand(4) = 100 * sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
            trials(3) = sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0);
            trials(4) = sum([Expt.Trials(:).bd] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
        else if strcmp(FileType, 'TWO')
                Grand(3) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir] ~= 0);
                Grand(4) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir] ~= 0);
                trials(3) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]<0 & [Expt.Trials(:).RespDir] ~= 0);
                trials(4) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd]>0 & [Expt.Trials(:).RespDir] ~= 0);
            else
                Grand(3) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0);
                Grand(4) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
                trials(3) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0);
                trials(4) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0);
                dxss(3) = mean([Expt.Trials([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]<0 & [Expt.Trials(:).RespDir] ~= 0).Id]);
                dxss(4) = mean([Expt.Trials([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id]>0 & [Expt.Trials(:).RespDir] ~= 0).Id]);
            end
        end
        g= []; d= []; t=[];
        for i = [ceil(size(dxValues,2)/2)+1:1:size(dxValues,2)-1];
            if strcmpi(FileType, 'BDpsych')
                g(i) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                d(i) = dxValues(i);
                t(i) = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            else if strcmp(FileType, 'TWO')
                    g(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                    d(i) = dxValues(i);
                    t(i) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                else
                    g(i) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                    d(i) = dxValues(i);
                    t(i) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                end
            end
        end
        Grand(5) = mean(g([ceil(size(dxValues,2)/2)+1:1:size(dxValues,2)-1]));
        dxss(5) = mean(d([ceil(size(dxValues,2)/2)+1:1:size(dxValues,2)-1]));
        trials(5) = sum(t([ceil(size(dxValues,2)/2)+1:1:size(dxValues,2)-1]));
        
        i=size(dxValues,2);
        if strcmpi(FileType, 'BDpsych')
            Grand(6) = 100 * sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            dxss(6) = dxValues(i);
            trials(6) = sum([Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
        else if strcmp(FileType, 'TWO')
                Grand(6) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                dxss(6) = dxValues(i);
                trials(6) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                Grand(12) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == -dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).bd] == -dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                if(NeuronNumber == 324 || NeuronNumber == 330)
                    Grand(12) = 100 * sum([Expt.Trials(:).dx] == dxValues(i-1) & [Expt.Trials(:).bd] == -dxValues(i-1) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i-1) & [Expt.Trials(:).bd] == -dxValues(i-1) & [Expt.Trials(:).RespDir] ~= 0);
                end
            else
                Grand(6) = 100 * sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] == ResponseToNegative) / sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
                dxss(6) = dxValues(i);
                trials(6) = sum([Expt.Trials(:).dx] == dxValues(i) & [Expt.Trials(:).RespDir] ~= 0);
            end
        end
        Grands{iN} = Grand;
        dxsz{iN} = dxss;
        Trialz{iN} = trials;
    else
        
    end
    
    if ShowIndividualPlots
        [slp,hp] = PsychSlop(NeuronName, StimulusType, FileType, 'dx', 1);
    end
    
    if SaveIndividualPlots
         print(hp, '-dpsc', '-r150', '-zbuffer', ['../figs/psych', '-', date, '-', filename, '.eps']);
    end
    
end

%% Grand 

figure(17653), 
i=0;
for ii = 1:length(Grands)
  if ~isempty(Grands{ii})
    if(Grands{ii}(1)<=Grands{ii}(6) || ii == 100)
        disp(strcat('PANIC!', num2str(ii)));
    else
        if(Grands{ii}(1) - Grands{ii}(6) >= 60) % performance above 80%
            if(~strcmp(FileType, 'TWO') || ( (Grands{ii}(3) - Grands{ii}(4)) >  mean([(Grands{ii}(1) - Grands{ii}(11)) , (Grands{ii}(12) - Grands{ii}(6))])))
                if(~strcmp(FileType, 'ABD') || (Grands{ii}(3) - Grands{ii}(4) > 0)) % ONLY bevahiorally EFFETIVE experiments 
                    disp(ii);
                    i = i+1;
                    g1(i) = Grands{ii}(1);
                    g2(i) = Grands{ii}(2);
                    g3(i) = Grands{ii}(3);
                    g4(i) = Grands{ii}(4);
                    g5(i) = Grands{ii}(5);
                    g6(i) = Grands{ii}(6);
                    if(g3(i)<=g4(i)), disp(strcat(num2str(ii) , ' - ' , num2str(g3(i)), ' - ' , num2str(g4(i)))), end
                    if strcmp(FileType, 'TWO')
                        g11(i) = Grands{ii}(11);
                        g12(i) = Grands{ii}(12);
                    end


                    id(i) = mean(abs(ids{ii}));

                    dx1(i) = dxsz{ii}(1);
                    dx2(i) = dxsz{ii}(2);
                    dx3(i) = dxsz{ii}(3);
                    dx4(i) = dxsz{ii}(4);
                    dx5(i) = dxsz{ii}(5);
                    dx6(i) = dxsz{ii}(6);

                    tr1(i) = Trialz{ii}(1);
                    tr2(i) = Trialz{ii}(2);
                    tr3(i) = Trialz{ii}(3);
                    tr4(i) = Trialz{ii}(4);
                    tr5(i) = Trialz{ii}(5);
                    tr6(i) = Trialz{ii}(6);
                end
            end
        end        
    end
  end
end
disp (i);

[mean(dx1), mean(dx2), mean(dx3), mean(dx4), mean(dx5), mean(dx6)]
[mean(tr1), mean(tr2), mean(tr3), mean(tr4), mean(tr5), mean(tr6)]
[sum(tr1), sum(tr2), sum(tr3), sum(tr4), sum(tr5), sum(tr6)]

if strcmp(FileType,'TWO')
    ggm = [mean(g11(~isnan(g11))), mean(g1), mean(g2(~isnan(g2))), mean(g3), mean(g5(~isnan(g5))), mean(g6), mean(g12(~isnan(g12))),];
    ggv = [std(g11(~isnan(g11))) /sqrt(i) ,std(g1) /sqrt(i) , std(g2(~isnan(g2))) /sqrt(i), std(g3) /sqrt(i), std(g5(~isnan(g5))) /sqrt(i), std(g6) /sqrt(i), std(g12(~isnan(g12))) /sqrt(i)];
else
    ggm = [mean(g1), mean(g2), mean(g3), mean(g5), mean(g6)];
    ggv = [std(g1) /sqrt(i) , std(g2) /sqrt(i), std(g3) /sqrt(i), std(g5) /sqrt(i), std(g6) /sqrt(i)];
end
errorbar(ggm, ggv, 'r');
hold on,
if strcmp(FileType,'TWO')
    ghm = [mean(g11(~isnan(g11))), mean(g1), mean(g2(~isnan(g2))), mean(g4), mean(g5(~isnan(g5))), mean(g6), mean(g12(~isnan(g12))),];
    ghv = [std(g11(~isnan(g11))) /sqrt(i) ,std(g1) /sqrt(i) , std(g2(~isnan(g2))) /sqrt(i), std(g4) /sqrt(i), std(g5(~isnan(g5))) /sqrt(i), std(g6) /sqrt(i), std(g12(~isnan(g12))) /sqrt(i)];
else
    ghm = [mean(g1), mean(g2), mean(g4), mean(g5), mean(g6)];
    ghv = [std(g1) /sqrt(i) , std(g2) /sqrt(i), std(g4) /sqrt(i), std(g5) /sqrt(i), std(g6) /sqrt(i)];
end
errorbar(ghm ,ghv, 'b');

%% Only correctly biased sessions
ggm = [mean(g1), mean(g2), mean(g3), mean(g5), mean(g6)];
ggv = [std(g1) /sqrt(i) , std(g2) /sqrt(i), std(g3) /sqrt(i), std(g5) /sqrt(i), std(g6) /sqrt(i)];
errorbar(ggm, ggv, 'r');
hold on,
ghm = [mean(g1), mean(g2), mean(g4), mean(g5), mean(g6)];
ghv = [std(g1) /sqrt(i) , std(g2) /sqrt(i), std(g4) /sqrt(i), std(g5) /sqrt(i), std(g6) /sqrt(i)];
errorbar(ghm ,ghv, 'b');




%%
load ../ScatterColors.mat
if strcmp(FileType,'TWO')
    figure, 
    scatter(((g1 - g11)), g3 - g4, [], colors(1:length(g3-g4),:), 'filled');    
    hold on, 
    scatter(((g12 - g6)), g3 - g4, [], colors(1:length(g3-g4),:), 'filled');
    scatter(mean([(g1 - g11);(g12 - g6)]), g3 - g4, [], colors(1:length(g3-g4),:), 'filled');
end

%% dx vs id

for i = 1: length(ids)
   idr(i) = range(ids{i});
   dxr(i) = range(dxs{i});
   idm(i) = mean(abs(ids{i}));
   dxm(i) = mean(abs(dxs{i}));
end


