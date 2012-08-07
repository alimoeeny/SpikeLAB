function fname = MakeFileName (MonkeyName, NeuronNumber, ClusterName, StimulusType, ExperimentType)
     if ClusterName(2) == 'e'
         fname = [MonkeyAb(MonkeyName), 'M', num2str(NeuronNumber, '%-04.3d'), '.', StimulusType,'.', ExperimentType, '.', ClusterName, '.mat'];
     else
         fname = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, StimulusType,'.', ExperimentType,'.mat');
     end