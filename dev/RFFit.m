function rf = RFFit(MonkeyName, NeuronNumber, ClusterName, ShowIndividualRFs)
% returns the rf as a vector [x, y of the center of receptive field, width (2*sd) of rf]
    DataPath = GetDataPath();
    % Do the ParaPos
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'PP');
%    filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.PP.mat');
    if (~strcmpi(ClusterName, '.c1.'))
        if (exist(filepath, 'file')~=2)
            ClusterName = '.c1.';
            %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.PP.mat');
            filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'PP');
        end
    end
    rfPP = [0 0 0];
    rfOP = [0 0 0];
    
    if (exist(filepath, 'file')==2)
        Neuron = load(filepath);
        Expt = Neuron.Expt;
       if isfield(Expt.Stimvals, 'Rx')
           rx = Expt.Stimvals.Rx;
           ry = Expt.Stimvals.Ry;
            for i = 1: length(Expt.Trials)
                Pps(i) = round(Expt.Trials(i).Pp*10.0)/10.0;
                helper(i,:) = [Pps(i), Expt.Trials(i).count];
            end
            values = unique(Pps);
            responses = [];
            for i = 1: length(values)
                responses(i) = mean(helper(helper(:,1) == values(i),2));
                distances(i) = values(i);
            end
        try
            cfun = fit(distances', responses', 'gauss1');
            rfPP = [cfun.b1 + rx, cfun.b1 + ry, cfun.c1 * 2];
        catch 
            rfPP = [0 0 0];
            disp('Fit Failed!')
        end
       else
           disp(['Something is wrong with this experiment!', filepath]);
        end
        
    else
        disp('Pp experiments not available!');
        rfPP = [0 0 0];
    end

    
    % Do the OrthoPos
    %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.OP.mat');
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'OP');
    if (~strcmpi(ClusterName, '.c1.'))
        if (exist(filepath, 'file')~=2)
            ClusterName = '.c1.';
            %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.OP.mat');
            filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'OP');
        end
    end
    
    if (exist(filepath, 'file')==2)
        Neuron = load(filepath);
        Expt = Neuron.Expt;
        if isfield(Expt.Stimvals, 'Rx')
           rx = Expt.Stimvals.Rx;
           ry = Expt.Stimvals.Ry;
            for i = 1: length(Expt.Trials)
                Ops(i) = round(Expt.Trials(i).Op*10.0)/10.0;
                helper(i,:) = [Ops(i), Expt.Trials(i).count];
            end
            values = unique(Ops);
            responses = [];
            for i = 1: length(values)
                responses(i) = mean(helper(helper(:,1) == values(i),2));
                distances(i) = values(i);
            end
        try
            cfun = fit(distances', responses', 'gauss1');
            rfOP = [cfun.b1 + rx, cfun.b1 + ry, cfun.c1 * 2];
        catch 
            rfPP = [0 0 0];
            disp('Fit Failed!')
        end
        else
            disp(['Something is wrong with this experiment!', filepath]);
        end
        
    else
        disp('Op experiments not available!');
        rfOP = [0 0 0];
    end
    
    % Do the Size
    %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.SZ.mat');
    filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'SZ');
    if (~strcmpi(ClusterName, '.c1.'))
        if (exist(filepath, 'file')~=2)
            ClusterName = '.c1.';
            %filename = strcat(MonkeyAb(MonkeyName), num2str(NeuronNumber, '%-04.3d'), ClusterName, 'rds.SZ.mat');
            filepath = MakeFilePath(MonkeyName, NeuronNumber, ClusterName, 'rds', 'SZ');
        end
    end
    
    if (exist(filepath, 'file')==2)
        Neuron = load(filepath);
        Expt = Neuron.Expt;
        values = sort(unique([Expt.Trials(:).sz]), 'descend');
        responses = [];
        for i = 1: length(values)
            responses(i) = mean([Expt.Trials([Expt.Trials(:).sz]==values(i)).count]);
        end
        SZ = values(responses == max(responses));
    else
        disp('SIZE experiments not available!');
        SZ = 0;
    end
    
    rf = [SZ rfOP rfPP];
    
    if  ShowIndividualRFs
        disp('Not implemented yet!')
    end