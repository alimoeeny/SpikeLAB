%clear;
%clc;

ShowSingleCells = 0; % 0 or 1


% % BDID
load('../AllBDIDNeuronsALL.mat');
AllNeurons = AllBDIDNeuronsALL;
clear AllBDIDNeuronsALL;





FileType = 'DT';
StimulusType = 'cylinder';



%Prep
DataPath = '/bgc/data/';
BinSize = 50;%50;
SmoothingBinSize = 1;%50
SmthKernel = gausswin(SmoothingBinSize);

filenamesforbruce = {};
TI=[];

%par
for iN= 1:length(AllNeurons) %[length(AllNeurons):-1:1], 1:length(AllNeurons)
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d')));
    
    DataPath = GetDataPath();
    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');

     if (exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename)))
        [eb, se, values, pD]  = PlotDT(MonkeyName, NeuronNumber, ClusterName, FileType, StimulusType, ShowSingleCells); 
        DT{iN} = {eb, se, values, pD};
     end
    
end


%%
for iN= 1:length(AllNeurons)
    if criteria(iN)
        if DT{iN}{4} == 1
            if DT{iN}{1}(end-1)~=0
                dt(iN,:) = DT{iN}{1};
                 v(iN,:) = DT{iN}{3};
            end
        else
            if DT{iN}{1}(end-1)~=0
                dt(iN,:) = DT{iN}{1}(length(DT{iN}{1}):-1:1);
                 v(iN,:) = DT{iN}{3};
            end
        end
    end
end


%%
figure, plot(mean(v(mean(dt,2)~=0,:)), mean(dt(mean(dt,2)~=0,:)));






