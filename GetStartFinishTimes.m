function [StartTime, FinishTime] = GetStartFinishTimes(filetype)

switch filetype
    case {'DID', 'DIDB'}
        StartTime = 500;
        FinishTime = 20000;
    case 'BDID'
        StartTime = 500;
        FinishTime = 5500;
    case 'ABD'
        StartTime = 500;
        FinishTime = 20000;
    case {'TWO'}
        StartTime = 500;
        FinishTime = 20000;
    otherwise 
        StartTime = 0;
        FinishTime = 0;

end

% StartTime = 5500; 
% if(isfield(Expt.Trials, 'dur'))
%     FinishTime = round(median([Expt.Trials(:).dur])) + 500;
% else
%     FinishTime = round(median([Expt.Trials(:).End] - [Expt.Trials(:).Start])) + 500;
% end
% if FinishTime>= StartTime
%     StartTime = 500;
% end
% if ~isempty(strfind(FileType, 'RID'))
%     FinishTime = 5000;
% end
