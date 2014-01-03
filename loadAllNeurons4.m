function [AllNeurons, FileType, StimulusType, StartTime, FinishTime] = loadAllNeurons4(ExptType)

switch ExptType
  case 'ABD'
    load('../AllABDNeurons.mat');
    AllNeurons = AllABDNeurons;
    FileType = 'ABD';
    StimulusType = 'cylinder';
    
  case 'DID'
    load('../AllDIDNeurons.mat');
    AllNeurons = AllDIDNeurons;
    FileType = 'DID';
    StimulusType = 'cylinder';
    StartTime  = 10000; %10000; % 6500; 
    FinishTime = 20000;
    
  case 'DTRW'
    load('../AllDTRWNeurons.mat');
    AllNeurons = AllDTRWNeurons;
    FileType = 'DTRW';
    StimulusType = 'cylinder';
    
  case 'TWO'
    load('../AllTWONeurons.mat');
    AllNeurons = AllTWONeurons;
    FileType = 'TWO';
    StimulusType = 'cylinder';
    StartTime  = 500; 
    FinishTime = 20000;
    
  case 'BDID'
    load('../AllBDIDNeuronsALL.mat');
    AllNeurons = AllBDIDNeuronsALL;
    FileType = 'BDID';
    StimulusType = 'cylinder';
    StartTime  = 10000; %10000; % 6500; 
    FinishTime = 20000;
    
  case 'DIDB'
    load('../AllDIDBNeurons.mat');
    AllNeurons = AllDIDBNeurons;
    FileType = 'DIDB';
    StimulusType = 'cylinder';
    
  case 'DRID'
    load('../AllDRIDNeurons.mat');
    FileType = 'DRID';
    AllNeurons = AllDRIDNeurons;
    StimulusType = 'rds';
    
  case 'SRID'
    load('../AllSRIDNeurons.mat');
    FileType = 'SRID';
    AllNeurons = AllSRIDNeurons;
    StimulusType = 'rds';
    
  case 'DPI'
    load('../AllPursuitNeurons.mat');
    AllNeurons = AllPursuitNeurons;
    FileType = 'DPI';
    StimulusType = 'cylinder';
    
  case 'DPIrds'
    load('../AllPursuitNeuronsrds.mat');
    AllNeurons = AllPursuitNeuronsrds;
    FileType = 'DPI';
    StimulusType = 'rds';

  case 'Psych'
    AllPsychData = importdata('AllDaedalusPsychDays.txt');
    AllNeurons = AllPsychData;
    clear AllPsychdata;
    FileType = 'psych';
    StimulusType = 'cylinder';

  case 'TWOBruce'
    load('../AllTWONeuronsFromBruce.mat');
    AllNeurons = AllNeurons;
    FileType = 'TWO';
    StimulusType = 'cylinder';
    StartTime  = 500; 
    FinishTime = 20000;
        
    
end