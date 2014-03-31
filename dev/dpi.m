
clear,
close all  
clc

cd /Users/moeenya/Dropbox/Projects/SpikeLAB/dev

DataPath = GetDataPath();
doitsquare = 0;

% DPI Cylinder
load('../AllPursuitNeurons.mat');
AllNeurons = AllPursuitNeurons;
clear AllPursuitNeurons;
FileType = 'DPI';
StimulusType = 'cylinder';
StartTime  = 500;
FinishTime = 20000;



for iN= [1:length(AllNeurons)] 
  NeuronNumber = AllNeurons(iN);
  [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronNumber); 
  disp(strcat('iN: ' ,num2str(iN) , ' , Neuron: ', num2str(NeuronNumber, '%-04.3d'), ' - ' , MonkeyName));
  TI(iN) = TuningIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, 'DT', [], doitsquare);
    
  filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', FileType,'.mat');
  Neuron = load(strcat(DataPath, MonkeyName, '/', num2str(NeuronNumber, '%-04.3d'), '/' ,filename));
  Expt = Neuron.Expt;
  fileNames{iN} = filename; 
    
  orS(iN) = ExperimentProperties(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType);
  [deltafxy, Speeds] = dpiDeltaFXY(Expt); 
 
  for sp = 1:length(Speeds)
    dxCm = sign([Expt.Trials(:).dx])==sign( TI(iN));
    dxCz = sign([Expt.Trials(:).dx])==0;
    dxCn = sign([Expt.Trials(:).dx])==sign(-TI(iN));
    if(isfield(Expt.Trials, 'dfx'))
        or = orS(iN);
        if ((orS(iN) > 180) || (orS(iN) < -180))
        or = - or;
        end
        purCp = sign([Expt.Trials(:).dfx])==sign(or);
        purCq = sign([Expt.Trials(:).dfx])==sign(-or);
    else
        or = orS(iN);
        if ((orS(iN) > 90) || (orS(iN) < -90))
            or = - or;
        end
        if(or==0), or = 1; end
        purCp = sign([Expt.Trials(:).dfy])==sign(-or);
      purCq = sign([Expt.Trials(:).dfy])==sign( or);      
    end
  
    ap = Expt.Trials(dxCm & purCp).count;
    az = Expt.Trials(dxCz & purCp).count;
    an = Expt.Trials(dxCn & purCp).count;
    bp = Expt.Trials(dxCm & purCq).count;
    bz = Expt.Trials(dxCz & purCq).count;
    bn = Expt.Trials(dxCn & purCq).count;
  end
  
  PIP(iN)    = (ap - bp) / (ap + bp);
  PIZero(iN) = (az - bz) / (az + bz);
  PINull(iN) = (an - bn) / (an + bn);
end