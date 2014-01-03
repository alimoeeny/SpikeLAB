function [E, cluster] = CutAndSave(DATA, varargin)    nosave = 0;    args = {};    j = 1;    while j <= length(varargin)        if strncmpi(varargin{j},'nosave',4)            nosave = 1;        else            args = {args{:} varargin{j}};        end        j = j+1;    end    if strcmp(DATA.autocutmode,'james')        [DATA, D]   = AllV.JamesAutoCut(DATA,'retrigger');        [AllVoltages, DATA] = AllV.BuildAllV(DATA, D.spk_inds(D.uids)', D.params.spk_pts);        outname = AllV.ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);        C = AllV.PlotTriggerHist(DATA,DATA.cluster);        if(0)            AllV.JamesAutoCut(DATA, 'reapply', C);        end        fprintf('Auto Cut %d events\n',max(DATA.uid));        DATA.savespikes = 1; %get the spikdi at least        DATA =  AllV.SaveClusters(DATA, outname,'quick');        cluster = C;        C.spkfile = AllV.SpkFileName(DATA,'auto');        E = C;        if ~exist(C.spkfile)            AllV.SaveSpikes(DATA,DATA.uid,C.spkfile);        end        DATA.savespikes = 0; %don't overwrite real spikes        set(DATA.toplevel,'UserData',DATA);        return;    end     [E, Scores, tmpdips, xy, details] = AllV.AutoCut(DATA, args{:});    DataClusters = AllV.mygetappdata(DATA,'Clusters');    if details.pctook > DATA.pcfit.took * 1.1        fprintf('PC fits %.2f vs %.2f\n',details.pctook, DATA.pcfit.took);    end    if ~isempty(Scores) && (~isfield(DATA,'TemplateScores') || E.plottype == 3)        DATA.TemplateScores = Scores;        DATA.tmpdips = tmpdips;    elseif E.newscores        DATA = get(DATA.toplevel,'UserData');    end    if E.angle ~= 0 %auto rotation caused by 1d/2d mismatch        DATA.xy{1} = xyrotate(xy(:,1),xy(:,2),E.angle);    else    DATA.xy{1} = xy;    end    DATA.ndxy = xy;    DATA.cboundary = E;    if AllV.IsTemplateCut(E)        DATA.plottype = 3;    else        DATA.plottype = E.plottype;    end    if DATA.watchplots%    AllV.PlotHistogram(DATA,E); %should be called in Classifyspikes    end    [cl, DATA.cluster, DATA.xy{1}] = AllV.ClassifySpikes(DATA,E);    DATA.cluster.auto = 1;        DATA.cluster.automode = DATA.autocutmode;    if  ~isfield(DATA.cluster,'dropi')        fprintf('no dropi calculated');    end    DATA.clid = cl.id;    DATA.nid = cl.nid;    DATA.clst = cl.clst;    DATA.MeanSpike = cl.MeanSpike';    DATA.cluster.MeanSpike = DATA.MeanSpike;    DATA.cluster.errs = DATA.errs;    DATA.cluster.starttime = DATA.cstarttime;    if isfield(DATA.Expt.Header,'ReadMethod')        DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;    else        DATA.cluster.exptreadmethod = 0;    end    % don't watndot do this unless saving, in which case its done below% This wipes out fields like savetimes%    DATA.Clusters{DATA.probe(1)} = DATA.cluster;            if nosave == 0        outname = AllV.ClusterFile(DATA.name,DATA.Expt,'auto','subdir',DATA.clustersubdir);        DATA =  AllV.SaveClusters(DATA, outname);        if length(DATA.savespkid) > length(DATA.xy{1})            fprintf('Save Spike Mismatch %d vs %d\n',length(DATA.savespkid),length(DATA.xy{1}));        end        if DATA.savespikes            AllV.SaveSpikes(DATA,DATA.savespkid,AllV.SpkFileName(DATA));        end    else            id = find(cl.clst(DATA.uid) > 1);% ClusterDetails records all event times, and the classification (clst)%Clusters just has the times of the classified  events = smallest file ffor%combine        DataClusters{DATA.probe(1)}.times = DATA.t(DATA.uid(id));    end    set(DATA.toplevel,'UserData',DATA);    cluster = DATA.cluster;        