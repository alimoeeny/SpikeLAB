% Eye track explore
clear, clc
DataPath = GetDataPath('server');

% % % % Eye Calib
% load('../AllEyeCallExpts.mat');
% AllNeurons = AllEyeCallExpts;
% clear AllEyeCallExpts;
% FileType = 'EyeCal';
% StimulusType = '';


% ABD
% load('/Users/ali/DropBox/Projects/BCode/AllABDNeurons.mat');
% AllNeurons = AllABDNeurons;
% clear AllABDNeurons;
% FileType = 'ABD';
% StimulusType = 'cylinder';
 
% % DID
[AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4('DID');
% load('../AllDIDNeurons.mat');
% AllNeurons = AllDIDNeurons;
% clear AllDIDNeurons;
% FileType = 'DID';
% StimulusType = 'cylinder';

% % TWO
% load('../AllTWONeurons.mat');
% AllNeurons = AllTWONeurons;
% clear AllTWONeurons;
% FileType = 'TWO';
% StimulusType = 'cylinder';

% % % DPI
% load('../AllPursuitNeurons.mat');
% AllNeurons = AllPursuitNeurons;
% clear AllPursuitNeurons;
% FileType = 'DPI';
% StimulusType = 'cylinder';

% % DPI rds
% load('../AllPursuitNeuronsrds.mat');
% AllNeurons = AllPursuitNeuronsrds;
% clear AllPursuitNeuronsrds;
% FileType = 'DPI';
% StimulusType = 'rds';


regenerateExsiting = 1;

%par
for iN= 1:length(AllNeurons), 
    Expt = [];
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , num2str(rem(now,1))));

    if (regenerateExsiting || ~(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'))==2))
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat'); 
        filename(findstr(filename, '..')) = '';
        
        filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);

        Neuron = load(filepath);
        Expt = Neuron.Expt;

%         if strcmpi(MonkeyName, 'icarus')
%             Expt = LoadEmData(Expt, 'lmonoc');
%             disp('Icarus - ');
%         else
%             Expt = LoadEmData(Expt, 'rmonoc');
%             disp('Daedalus - ');
%         end

        % Just load all the Eye movement data
        Expt = LoadEmData(Expt)
        
        if(isfield(Expt.Trials, 'EyeData'))
            %r = myParSave(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'), Expt); 
            filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPath);
            filepath = [filepath(1:end-3), 'EYE.mat'];
            r = myParSave(filepath, Expt);
            if r ~= 0, 
                disp(' W A I T ! something is not quite right here!'); 
            end
        end
    end
end