function testAGenerateTheMockData

maxTrialSpikeCount = 100;
trialLength = 500;
start = 0;


%% rds.DT -> TI 

% %%%%%%% ALL ZERO 
clear Expt;
Expt.Header.Name = 'test001.c1.rds.DT.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:73
    Expt.Trials(i).count = 0;
    Expt.Trials(i).dx = 0;
    Expt.Trials(i).Start = 0;
    Expt.Trials(i).End = 0;
    Expt.Trials(i).Spikes = [];
end

pth = '/bgc/data/test/001/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI >> 0
clear Expt;
Expt.Header.Name = 'test002.c1.rds.DT.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:97
    r = randi(2);
    switch r
        case 1
          count = randi(maxTrialSpikeCount/10);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/002/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI << 0
clear Expt;
Expt.Header.Name = 'test003.c1.rds.DT.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:97
    r = randi(2);
    switch r
        case 1
          count = randi(maxTrialSpikeCount);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount/10);
          dx = -0.01;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/003/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI ~ 0
clear Expt;
Expt.Header.Name = 'test004.c1.rds.DT.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:157
    r = randi(2);
    switch r
        case 1
          count = randi(maxTrialSpikeCount);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/004/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');

% %%%%%%% GENERIC TI > 0 with zero dx
clear Expt;
Expt.Header.Name = 'test005.c1.rds.DT.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:157
    r = randi(3);
    switch r
        case 1
          count = randi(maxTrialSpikeCount/10);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
        case 3
          count = randi(maxTrialSpikeCount/5);
          dx = 0.0;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/005/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');



%% rds.DPI -> PI 

% %%%%%%% ALL ZERO 
clear Expt;
Expt.Header.Name = 'test001.c1.rds.DPI.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:49
    Expt.Trials(i).count = 0;
    Expt.Trials(i).dx = 0;
    Expt.Trials(i).Start = 0;
    Expt.Trials(i).End = 0;
    Expt.Trials(i).Spikes = [];
    Expt.Trials(i).dfy = 3;
    Expt.Trials(i).fy = 0;
end

pth = '/bgc/data/test/001/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI >> 0
clear Expt;
Expt.Header.Name = 'test002.c1.rds.DPI.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:97
    r = randi(12);
    switch r
        case 1
          count = randi(maxTrialSpikeCount/10);
          dx = 0.01;
          dfy = 3;
          fy = 0;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.0;
          dfy = 3;
          fy = 0;
        case 3
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
          dfy = 3;
          fy = 0;
        case 4
          count = randi(maxTrialSpikeCount/5);
          dx = 0.01;
          dfy = -3;
          fy = 0;
        case 5
          count = randi(maxTrialSpikeCount);
          dx = -0.0;
          dfy = -3;
          fy = 0;
        case 6
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
          dfy = -3;
          fy = 0;
        case 7
          count = randi(maxTrialSpikeCount/10);
          dx = 0.01;
          dfy = 2;
          fy = 0;
        case 8
          count = randi(maxTrialSpikeCount);
          dx = -0.0;
          dfy = 2;
          fy = 0;
        case 9
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
          dfy = 2;
          fy = 0;
        case 10
          count = randi(maxTrialSpikeCount/5);
          dx = 0.01;
          dfy = -2;
          fy = 0;
        case 11
          count = randi(maxTrialSpikeCount);
          dx = -0.0;
          dfy = -2;
          fy = 0;
        case 12
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
          dfy = -2;
          fy = 0;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
    Expt.Trials(i).dfy = dfy;
    Expt.Trials(i).fy = fy;
end

pth = '/bgc/data/test/002/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI << 0
clear Expt;
Expt.Header.Name = 'test003.c1.rds.DPI.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:97
    r = randi(2);
    switch r
        case 1
          count = randi(maxTrialSpikeCount);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount/10);
          dx = -0.01;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/003/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');


% %%%%%%% GENERIC TI ~ 0
clear Expt;
Expt.Header.Name = 'test004.c1.rds.DPI.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:157
    r = randi(2);
    switch r
        case 1
          count = randi(maxTrialSpikeCount);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/004/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');

% %%%%%%% GENERIC TI > 0 with zero dx
clear Expt;
Expt.Header.Name = 'test005.c1.rds.DPI.mat';
Expt.Stimvals.or = 0;
Expt.Stimvals.et = 'dx';

for i = 1:157
    r = randi(3);
    switch r
        case 1
          count = randi(maxTrialSpikeCount/10);
          dx = 0.01;
        case 2
          count = randi(maxTrialSpikeCount);
          dx = -0.01;
        case 3
          count = randi(maxTrialSpikeCount/5);
          dx = 0.0;
    end
    
    spikes = [];            
    for s = 1:count,  spikes(s) = 10 * s * trialLength / count + randi(100); end
    spikes = sort(spikes);

    
    Expt.Trials(i).count = count;
    Expt.Trials(i).dx = dx;
    Expt.Trials(i).Start = start;
    Expt.Trials(i).End = start + trialLength;
    Expt.Trials(i).Spikes = spikes;
end

pth = '/bgc/data/test/005/';
if ~(exist(pth,'dir')==7), mkdir(pth); end
save(strcat(pth,Expt.Header.Name), 'Expt');



