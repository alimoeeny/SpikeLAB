function [s] = MonkeyAb(MonkeyName)

mn = MonkeyName;

if strcmpi(mn, 'ICARUS')
    s = 'ic';
else if strcmpi(MonkeyName, 'DAE')
    s = 'dae';
else if strcmpi(MonkeyName, 'TEST')
    s = 'test';
    end
    end
end