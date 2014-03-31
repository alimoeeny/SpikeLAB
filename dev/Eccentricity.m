function ex = Eccentricity(NeuronName, StimulusType, FileType, DataPaths)
   [MonkeyName, NeuronNumber, ClusterName] = NeurClus(NeuronName);
   filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, StimulusType, FileType, DataPaths);
   Neuron = load(filepath);
   Expt = Neuron.Expt;
   
   xo = Expt.Stimvals.xo;
   yo = Expt.Stimvals.yo;
   if(isempty(xo))
       xo = mean([Expt.Trials.xo]);
   end
   if(isempty(yo))
       yo = mean([Expt.Trials.yo]);
   end
   
   ex = sqrt(xo * xo + yo * yo);