function [ results ] = calibrateICC(oldICCFilePath, newICCFilePath, intendedLuminances, measuredLuminances)
% reads the current (old) lookup table from the old ICC file
% fits an offseted gamma function to the measurements like
% y = min(M) + (range(M) * x ^ gamma) WHERE x is intended_luminances /
% max(intended luminances) and M is the measured_luminance and y is the
% measured_luminances
% Ali 

results = 'started';
assert(sum(size(measuredLuminances)==size(intendedLuminances))==2);

targetKernel = linspace(0,256,256);

a = min(measuredLuminances);
b = range(measuredLuminances);
fitdescription = fittype([num2str(a) ' + ' num2str(b) '  * (x^c)']);
disp([fitdescription]);

%[cfun,gof,output] = fit(intendedLuminances', measuredLuminances', fitdescription);

[fo,gof] = fit(intendedLuminances', measuredLuminances', 'poly3');

%measuredGamma = cfun.c;

%measurmentSimulation = a + b * linspace(0,1,256) .^measuredGamma;

x = linspace(0,1,256);
measurmentSimulation = fo.p1*x.^3 + fo.p2*x.^2 + fo.p3.*x + fo.p4;

figure(1182), clf, hold on, 
plot(linspace(0,1,256).*256, measurmentSimulation);
plot(intendedLuminances*256, measuredLuminances, 'r');

mycc = iccread(oldICCFilePath);

oldLookupTable = [];
for i = 1:size(mycc.PrivateTags,1)
    if(strcmp(mycc.PrivateTags{i,1}, 'vcgt'))
        oldLookupTable = mycc.PrivateTags{i,2}(19:2:530);
        ptd = mycc.PrivateTags{i,2};
        vgctIndex = i;
        break
    end
end
if(isempty(oldLookupTable))
    disp('Cant read the lookuptable from ICC file');
    return;
end

measurmentSimulationScaled  = measurmentSimulation ./ max(measurmentSimulation) .* max(targetKernel); ;
y = double(oldLookupTable) .* targetKernel ./ measurmentSimulationScaled;

figure(1919), clf, hold on
plot(double(oldLookupTable), 'b');
plot(targetKernel, 'r');
plot(measurmentSimulationScaled, 'g');
plot(y, 'm');

n16y = double(y);
n16y = (n16y - min(n16y)) ./ range(n16y);
n16y = n16y * 2^16;
n16y = uint16(n16y);
n16y = sixteen2eight(n16y);

% 1:18 are the strings don't change them

% 19:530 are the fisrt gun shoud be RED
%ptd(19:2:530) = y; %high bits
%ptd(20:2:530) = 0; %low bits
ptd(19:530) = n16y; 

% 531:1042 are the second gun shoud be GREEB
%ptd(531:2:1042) = y; %high bits
%ptd(532:2:1042) = 0; %low bits
ptd(531:1042) = n16y;

% 1043:1554 are the third gun shoud be BLUE
%ptd(1043:2:1554) = y; %high bits
%ptd(1044:2:1554) = 0; %low bits
ptd(1043:1554) = n16y;

mycc.PrivateTags{vgctIndex,2} = ptd;

% the name that appears in the preferences window in the list of profiles
%mycc.Description.String = ['ColorLCDGamma' num2str(measuredGamma)  '.', num2str(now)];
newName = ['Ali' char('A' + randi(24,30,1))'];
mycc.Description.String = newName;

rc = 1;
for ri = 1:2:35
    if(ri<=length(newName))
        mycc.PrivateTags{1,2}(37+ri) = uint8(newName(rc));
    else
        mycc.PrivateTags{1,2}(37+ri) = uint8(0);
    end
    rc = rc + 1;
end

iccwrite(mycc, [newICCFilePath newName '.icc']);
results = 'success';

end

