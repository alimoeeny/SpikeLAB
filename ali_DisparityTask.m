function err = ali_DisparityTask(varargin)

arguments = parsevarargin(varargin);

basedir = '/local/AliManStims/';
if arguments.isKey('basedir')
    basedir = arguments('basedir');
end
% 
% [ST,~] = dbstack;
% experimentname = ST.name;

numberofrepetitions = 4;

%st
stimulustype = 'rds';
%et
experimentalvariable = 'dx';
%nt
numberofregularconditions = 13;
%mt
meanstimvalue = 0;
%in
increment = 0.2;


%DEFAULT VALUES
xo =-10;
yo =-10;
sz = 12;
jv = 10;
nf = 240;


baselinesetup = containers.Map;
baselinesetup('fs') = 0.3;
baselinesetup('st') = stimulustype;
baselinesetup('fw') = 1.3;
baselinesetup('xo') = xo;
baselinesetup('yo') = yo;
baselinesetup('sz') = sz;
baselinesetup('jv') = jv;
baselinesetup('Rx') = xo;
baselinesetup('Ry') = yo;
baselinesetup('nf') = nf;
baselinesetup('sl') = 0;
baselinesetup('co') = 1;
baselinesetup('ce') = 1;
baselinesetup('op') = '+exm+sq';
baselinesetup('!mat') = [experimentname];

filename = [basedir experimentname '.stm'];
err = ali_savemap2file(baselinesetup, filename);
if err ~= ''
    disp('SOMETHING WHENT WRONG HERE SAVEING BASELINE EXPERIMENT')
    return
end

if exist([basedir experimentname],'dir')
    delete([basedir experimentname '/*']);
    rmdir([basedir experimentname]);
end
[SUCCESS,~,~] = mkdir([basedir experimentname]);
if SUCCESS == 0
    disp('ERROR CREATING THE stim DIRECTORY, CANNOT CONTINUE');
    err = 'ERROR CREATING THE stim DIRECTORY';
    return
end


stimcounter = 0; % this is the suffix of the stim files and is zero based

stimvalues = -floor(numberofregularconditions / 2) * increment + meanstimvalue :increment: floor(numberofregularconditions / 2) * increment + meanstimvalue;
for i = stimvalues
    stim = containers.Map;
    stim(experimentalvariable) = i;
    
    filename = [basedir experimentname '/stim' num2str(stimcounter)];
    stimcounter = stimcounter + 1;
    err = ali_savemap2file(stim, filename);
    if err ~= ''
        disp('SOMETHING WHENT WRONG HERE SAVEING STIM FILE')
        return
    end
end


oddballstim = containers.Map;
oddballstim('ce') = 0;
oddballstim('st') = 'rds';

filename = [basedir experimentname '/stim' num2str(stimcounter)];
stimcounter = stimcounter + 1;
err = ali_savemap2file(oddballstim, filename);
if err ~= ''
    disp('SOMETHING WHENT WRONG HERE SAVEING BASELINE EXPERIMENT')
    return
end

% - - - - - - - - - - - - STIM ORDER file
stims = 0:stimcounter-1;
% to make the oddball apear twice as much as others
stims(end+1) = stims(end);

overallstimorder = [];
for rpt = 1:numberofrepetitions
    overallstimorder = [overallstimorder stims(randperm(length(stims)))];
end
filename = [basedir experimentname '/stimorder'];
dlmwrite(filename, overallstimorder, ' ');

disp(['Total stim files created:', num2str(stimcounter)]);
err = '';
