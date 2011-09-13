function [rv] = doRDSRLS(iv) 
% disparity tuning in MT for rls vs rds

DataPath = GetDataPath();

load ../AllrdsDXST.mat
AllNeurons = AllrdsDXST;
clear AllrdsDXST;
StimulusType = 'rds';
FileType = 'DXST';
StartTime = 500;
FinishTime = 5000;

NeuronsRange = 1:length(AllNeurons);

for iN = NeuronsRange
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(AllNeurons(iN)); 
    disp([num2str(iN, '%-3.3d'), ' - Neuron: ', num2str(NeuronNumber, '%-04.3d')]);

    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
    Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
    Expt = Neuron.Expt;
    fileNames{iN} = filename;
    
    stS = expMining(Expt, [], 'st');
    dxS = expMining(Expt, [], 'dx');
    
%     FRs(1, :) = expFR(Exp, stS);
%     FRs(2, :) = expFR(Exp, dxS);
    R{iN} = expFR(Expt, StartTime, FinishTime, stS, dxS);
    
    %ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, 'bo');
    h = figure(1235+iN), errorbar([R{iN}{4}', R{iN}{4}'], R{iN}{1}', R{iN}{2}')
    set(h, 'LineWidth', 2);
    set(cga, 'XGrid', 'on');
    title('Disparity tuning for rds vs rls');
    legend('rds', 'rls')
end

rv = R;
