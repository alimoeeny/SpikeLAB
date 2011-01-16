function [s] = MonkeyAb(MonkeyName)

if strcmpi(MonkeyName, 'icarus')
    s = 'ic';
else if strcmpi(MonkeyName, 'dae')
    s = 'dae';
    end
end