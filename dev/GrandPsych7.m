clc;
clear;
DataPath = GetDataPath();


ShowIndividualPlots = 0; % 0 or 1
SaveIndividualPlots = 0; % 0 or 1
[AllNeurons, FileType, StimulusType] = loadAllNeurons4('TWOBruce');

%load ~/Downloads/twoall.mat

% THE REAL THING
%refrenceDxs = [-0.0335   -0.0067   -0.0033         0    0.0033    0.0067    0.0335];
% THE REAL THING  wiht abs(TI)>0.05
if strcmp(FileType, 'TWO')
    %refrenceDxs = [-0.0126   -0.0030   -0.0023   -0.0001    0.0025    0.0034    0.0126];
    % Bruces Selection
    refrenceDxs = [-0.0122   -0.003   -0.002   0.000    0.002    0.0030    0.0122];
    %refrenceDxs = [-0.0122   -0.0025   -0.0019   -0.0001    0.0022    0.0030    0.0122];
else
    % abs(TI)>0.1 refrenceDxs = [-0.0280   -0.0051   -0.0021    0.00    0.0031    0.0098    0.0287];
    refrenceDxs = [-0.0184   -0.0035   -0.0017    0.0000    0.0018    0.0043    0.0185];
end


vivals = zeros(length(AllNeurons),7);
for iN= 1:length(AllNeurons), %iN= [20: 25] %length(AllNeurons)] %[1:20] % 1:length(AllNeurons),
    NeuronName = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronName);
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));
    
    filename = MakeFileName(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
    
    TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, []);
     
        Neuron = load(filepath);
        Expt = Neuron.Expt;
        dxValues = unique([Expt.Trials(:).dx]);
        
        if(mean([Expt.Trials([Expt.Trials(:).dx]>0).RespDir])>0)
            ResponseToPositive = 1;
            ResponseToNegative = -1;
        else
            ResponseToPositive = -1;
            ResponseToNegative = 1;
        end
        
        Grand = zeros(7,1);
        dxss = [];
        trials = zeros(7,1);
        
        Grands(iN,1:7) = 0;
        GrandsId(iN, 1:2) = 0;
        GrandsFlip(iN, 1:2) = 0;
        
        valuesLengths(iN) = length(dxValues);
        
        if (length(dxValues)>14)
            continue
        end
        
        
        
        if ((dxValues(2)<refrenceDxs(1)) || iseven(length(dxValues)))
            disp('we are not gonna use this one!');
        else
            for dxi = 1:length(dxValues)
                if(dxValues(dxi)~=0)
                    if (length(dxValues)==7)
                        dxii = dxi;
                    else if (length(dxValues)==5)
                            switch dxi
                                case 1, dxii = 1;
                                case 2, dxii = 2;
                                case 3, dxii = 4;
                                case 4, dxii = 6;
                                case 5, dxii = 7;
                            end
                        else if (length(dxValues)==9)
                                switch dxi
                                    case 1, dxii = 1;
                                    case 2, dxii = 2;
                                    case 3, dxii = 2;
                                    case 4, dxii = 3;
                                    case 5, dxii = 4;
                                    case 6, dxii = 5;
                                    case 7, dxii = 6;
                                    case 8, dxii = 6;
                                    case 9, dxii = 7;
                                end
                            else if (length(dxValues)==11)
                                    switch dxi
                                        case 1, dxii = 1;
                                        case 2, dxii = 2;
                                        case 3, dxii = 2;
                                        case 4, dxii = 3;
                                        case 5, dxii = 3;
                                        case 6, dxii = 4;
                                        case 7, dxii = 5;
                                        case 8, dxii = 5;
                                        case 9, dxii = 6;
                                        case 10,dxii = 6;
                                        case 11,dxii = 7;
                                    end
                                else if (length(dxValues)==13)
                                        switch dxi
                                            case 1, dxii = 1;
                                            case 2, dxii = 2;
                                            case 3, dxii = 2;
                                            case 4, dxii = 3;
                                            case 5, dxii = 3;
                                            case 6, dxii = 3;
                                            case 7, dxii = 4;
                                            case 8, dxii = 5;
                                            case 9, dxii = 5;
                                            case 10,dxii = 5;
                                            case 11,dxii = 6;
                                            case 12,dxii = 6;
                                            case 13,dxii = 7;
                                        end
                                    else
                                        if (dxValues(dxi)<0)
                                            dxii = find(abs(refrenceDxs-dxValues(dxi))==min(abs(refrenceDxs-dxValues(dxi))),1, 'first');
                                        else
                                            dxii = find(abs(refrenceDxs-dxValues(dxi))==min(abs(refrenceDxs-dxValues(dxi))),1, 'last');
                                        end
                                        
                                        if(refrenceDxs(dxii)==0 && dxValues(dxi)~=0)
                                            dxii = dxii + (1 * sign(dxValues(dxi)));
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    
                    
                    if(vivals(iN, dxii))
                        vivals(iN, dxii) = mean([vivals(iN, dxii), dxValues(dxi)]);
                    else
                        vivals(iN, dxii) = dxValues(dxi);
                    end
                    
                    if((dxii == 4) || (dxii == 8))
                        debug =1 ;
                    end
                    
                    if((dxii>4 && dxValues(dxi)<0) || (dxii<4 && dxValues(dxi)>0))
                        debug = 1;
                    end
                    
                    
                    if(Grand(dxii))
                        Grand(dxii)  = mean([Grand(dxii), ...
                            100 * sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] ~= 0)]);
                    else
                        Grand(dxii)  = 100 * sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] ~= 0);
                    end
                    
                    
                    if(trials(dxii))
                        trials(dxii) =  mean([trials(dxii), ...
                            sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] ~= 0)]);
                    else
                        trials(dxii) =  sum([Expt.Trials(:).dx] == dxValues(dxi) & [Expt.Trials(:).RespDir] ~= 0);
                    end
                else
                    
                   
                end
                
            end
            
            if (length(dxValues)==5)
                vivals(iN, 3) = vivals(iN, 2);
                vivals(iN, 5) = vivals(iN, 6);
                Grand(3) = Grand(2);
                Grand(5) = Grand(6);
                trials(3) = trials(2);
                trials(5) = trials(6);
            end
            
            if strcmp(FileType, 'TWO')
                GrandsId(iN, 1) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] > 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] > 0 & [Expt.Trials(:).RespDir] ~= 0);
                GrandsId(iN, 2) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] < 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] < 0 & [Expt.Trials(:).RespDir] ~= 0);
                GrandsFlip(iN, 1) = 100 * sum([Expt.Trials(:).dx] == (-1 * [Expt.Trials(:).bd]) & [Expt.Trials(:).bd] > 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] > 0 & [Expt.Trials(:).RespDir] ~= 0);
                GrandsFlip(iN, 2) = 100 * sum([Expt.Trials(:).dx] == (-1 * [Expt.Trials(:).bd]) & [Expt.Trials(:).bd] < 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] < 0 & [Expt.Trials(:).RespDir] ~= 0);
                TrialsId(iN, 1) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] > 0 & [Expt.Trials(:).RespDir] ~= 0);
                TrialsId(iN, 2) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).bd] < 0 & [Expt.Trials(:).RespDir] ~= 0);
            else
                GrandsId(iN, 1) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] > 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] > 0 & [Expt.Trials(:).RespDir] ~= 0);
                GrandsId(iN, 2) = 100 * sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] < 0 & [Expt.Trials(:).RespDir] == ResponseToPositive) / sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] < 0 & [Expt.Trials(:).RespDir] ~= 0);
                TrialsId(iN, 1) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] > 0 & [Expt.Trials(:).RespDir] ~= 0);
                TrialsId(iN, 2) = sum([Expt.Trials(:).dx] == 0 & [Expt.Trials(:).Id] < 0 & [Expt.Trials(:).RespDir] ~= 0);
            end
            
            if strcmp(FileType, 'TWO')
                ids(iN,1:2) = [-1 * mean(abs(unique([Expt.Trials(:).bd]))) , mean(abs(unique([Expt.Trials(:).bd])))];
            else    
                ids(iN,1:2) = [-1 * mean(abs(unique([Expt.Trials(:).Id]))) , mean(abs(unique([Expt.Trials(:).Id])))]; %unique([Expt.Trials(:).Id]);
            end
            %Grand(isnan(Grand)) = 0;
            Grands(iN,1:7) = Grand;
            Trialz(iN,1:7) = trials;
        end
        
end


%% Grand

%%
figure(18854), clf, hold on,

xxm = mean(vivals);
ggm = mean(Grands(mean(Grands')~=0, :));
ggv = std(Grands(mean(Grands')~=0, :)) ./ sqrt(length(AllNeurons));
errorbar(xxm, ggm, ggv, 'r');

ggmId = mean(GrandsId(mean(GrandsId')~=0, :));
ggvId = std(GrandsId(mean(GrandsId')~=0, :)) ./ sqrt(length(AllNeurons)); % sqrt(TrialsId(mean(Grands')~=0,:));
errorbar([0.001 -0.001], ggmId, ggvId, 'b');

%%
figure(12256), clf, hold on,

xxm = mean(vivals(abs(TI)>0.05, :));
ggm = mean(Grands(mean(Grands')~=0& abs(TI)>0.05, :),1);
ggv = std(Grands(mean(Grands')~=0 & abs(TI)>0.05, :),[],1) ./ sqrt(sum(abs(TI)>0.05));
errorbar(xxm, ggm, ggv, 'r');

ggmId = mean(GrandsId(mean(GrandsId')~=0 & abs(TI)>0.05, :));
ggvId = std(GrandsId(mean(GrandsId')~=0 & abs(TI)>0.05, :)) ./ sqrt(sum(abs(TI)>0.05)); % sqrt(TrialsId(mean(Grands')~=0,:));
errorbar([0.001 -0.001], ggmId, ggvId, 'b');


%% Bruce's selection
load ~/Dropbox/Projects/SpikeLAB/AllTWONeuronsFromBruce.mat BrucesSelection
figure(42256), clf, hold on,

xxm = mean(vivals(BrucesSelection, :));
ggm = mean(Grands(mean(Grands')~=0& BrucesSelection, :));
ggv = std(Grands(mean(Grands')~=0 & BrucesSelection, :)) ./ sqrt(sum(BrucesSelection));
errorbar(xxm, ggm, ggv, 'r');

ggmId = mean(GrandsId(mean(GrandsId')~=0 & BrucesSelection, :));
ggvId = std(GrandsId(mean(GrandsId')~=0 & BrucesSelection, :)) ./ sqrt(sum(BrucesSelection)); % sqrt(TrialsId(mean(Grands')~=0,:));
errorbar([0.001 -0.001], ggmId, ggvId, 'b');

ggmFlip = mean(GrandsFlip(mean(GrandsFlip')~=0 & BrucesSelection, :));
ggvFlip = std(GrandsFlip(mean(GrandsFlip')~=0 & BrucesSelection, :)) ./ sqrt(sum(BrucesSelection)); % sqrt(TrialsId(mean(Grands')~=0,:));
errorbar([0.02 -0.02], ggmFlip, ggvFlip, 'g');

%% for cftool

x = xxm([1 2 3 5 6 7]);
 g = ggm([1 2 3 5 6 7]);
 xf = x;
 gf = g / 100;
 
 cftool

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is the fit
fit = 0.5 + erf((x - a)/(b * sqrt(2)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
