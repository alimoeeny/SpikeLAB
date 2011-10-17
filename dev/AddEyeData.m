% Eye track explore
clear, clc
DataPath = GetDataPath();

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
load('../AllPursuitNeurons.mat');
AllNeurons = AllPursuitNeurons;
clear AllPursuitNeurons;
FileType = 'DPI';
StimulusType = 'cylinder';

%par
parfor iN= 1:length(AllNeurons), 
    Expt = [];
    NeuronNumber = AllNeurons(iN);
    [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
    disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , num2str(rem(now,1))));

    if ~(exist(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' , MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'))==2)
        filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat'); 
        filename(findstr(filename, '..')) = '';
        Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
        Expt = Neuron.Expt;

        if strcmpi(MonkeyName, 'icarus')
            Expt = LoadEmData(Expt, 'lmonoc');
            disp('Icarus - ');
        else
            Expt = LoadEmData(Expt, 'rmonoc');
            disp('Daedalus - ');
        end
        if(isfield(Expt.Trials, 'EyeData'))
            r = myParSave(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/', MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.Eye.mat'), Expt); 
            if r ~= 0, disp(' W A I T ! something is not quite right here!'); end
        end
    end
end