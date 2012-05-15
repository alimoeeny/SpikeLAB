clear


filename = '/bgc/data/dae/062/dae062.c1.cylinder.DID.mat';
disp(filename);
load(filename);


filename = [filename(1:end-3), 'txt'];
dlmwrite(filename, 'ConvertDate', 'delimiter', '');
dlmwrite(filename, datestr(date), 'delimiter', '', '-append');

fns = fieldnames(Expt.Header);
for i = 1: length(fns);
    %disp(Expt.Header.(fns{i}));
    if (sum(strcmpi(fns{i},{'Clusters', 'Peninfo', 'SpkStats', 'Cluster'}))==0)
        dlmwrite(filename, fns{i}, 'delimiter', '', '-append');
        S = [];
        if (isempty(Expt.Header.(fns{i})))
            dlmwrite(filename, ' ', '-append');
        else
            if (iscell(Expt.Header.(fns{i})))
                if (isstr(Expt.Header.(fns{i}){:}))
                    S = strrep(Expt.Header.(fns{i}), char(10),'');
                end
            else if (isstr(Expt.Header.(fns{i})(:)))
                        S = strrep(Expt.Header.(fns{i}), char(10),'');
                 end
            end
            if ~isempty(S)
                dlmwrite(filename, S, 'delimiter', '', '-append');
            else
                dlmwrite(filename, Expt.Header.(fns{i}), 'precision', 15, '-append');
            end
        end
    end
end

fns = fieldnames(Expt.Stimvals);
for i = 1: length(fns);
    %disp(Expt.Stimvals.(fns{i}));
    if (sum(strcmpi(fns{i},{}))==0) %'Clusters', 'Peninfo', 'SpkStats', 'Cluster'}))==0)
        dlmwrite(filename, fns{i}, 'delimiter', '', '-append');
        S = [];
        if (isempty(Expt.Stimvals.(fns{i})))
            dlmwrite(filename, ' ', '-append');
        else
            if (iscell(Expt.Stimvals.(fns{i})))
                if (isstr(Expt.Stimvals.(fns{i}){:}))
                    S = strrep(Expt.Stimvals.(fns{i}), char(10),'');
                end
            else if (isstr(Expt.Stimvals.(fns{i})(:)))
                        S = strrep(Expt.Stimvals.(fns{i}), char(10),'');
                 end
            end
            if ~isempty(S)
                dlmwrite(filename, S, 'delimiter', '', '-append');
            else
                dlmwrite(filename, Expt.Stimvals.(fns{i}), 'precision', 15, '-append');
            end
        end
    end
end



fns = fieldnames(Expt.Trials);
for t = 1: length(Expt.Trials);
%    dlmwrite(filename, 'Trial=>', 'delimiter', '', '-append');
%    dlmwrite(filename, t, '-append');
    for i = 1: length(fns);
        %disp(fns{i});
        %disp(Expt.Trials(t).(fns{i}));
        if (sum(strcmpi(fns{i},{'Events'}))==0) %'Clusters', 'Peninfo', 'SpkStats', 'Cluster'}))==0)
            dlmwrite(filename, ['Trial=>', num2str(t, '%-04.4d'), fns{i}], 'delimiter', '', '-append');
            S = [];
            if (isempty(Expt.Trials(t).(fns{i})))
                dlmwrite(filename, ' ', '-append');
            else
                if (iscell(Expt.Trials(t).(fns{i})))
                    if (isstr(Expt.Trials(t).(fns{i}){:}))
                        S = strrep(Expt.Trials(t).(fns{i}), char(10),'');
                    end
                else if (isstr(Expt.Trials(t).(fns{i})(:)))
                        S = strrep(Expt.Trials(t).(fns{i}), char(10),'');
                    end
                end
                if ~isempty(S)
                    dlmwrite(filename, S, 'delimiter', '', '-append');
                else
                    if (sum(strcmpi(fns{i},{'Spikes', 'OSpikes', 'Ocodes'}))==0) 
                        dlmwrite(filename, Expt.Trials(t).(fns{i}), 'precision', 15, '-append');
                    else
                        dlmwrite(filename, Expt.Trials(t).(fns{i})', 'precision', 15, '-append');
                    end
                end
            end
        end
    end
end






