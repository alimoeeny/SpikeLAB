function [ results ] = calibrateICCStereo(oldICCFilePath, outputPath, intendedLuminances, measuredLuminancesRightBlue, measuredLuminancesLeftRed)
% oldICCFilePath the path to the icc profile used to do the measurements
% outputpath is the path to the directory where the new icc file will be written
% intendedLuminances is usually somehting like: linspace(0,1,21) for a 21level graylevels measurment
% measuredLuminancesRightBlue,  measuredLuminancesLeftRed are the actual measured values (photometer) for right and left monitors/screens
% reads the current (old) lookup table from the old ICC file
% fits an offseted gamma function to the measurements like
% y = min(M) + (range(M) * x ^ gamma) WHERE x is intended_luminances /
% max(intended luminances) and M is the measured_luminance and y is the
% measured_luminances
% Ali 

results = 'started';
assert(sum(size(measuredLuminancesRightBlue)==size(intendedLuminances))==2);
assert(sum(size(measuredLuminancesLeftRed)==size(intendedLuminances))==2);

targetKernel = linspace(0,256,256);

aRight = min(measuredLuminancesRightBlue);
bRight = range(measuredLuminancesRightBlue);
aLeft = min(measuredLuminancesLeftRed);
bLeft = range(measuredLuminancesLeftRed);
fitdescriptionRight = fittype([num2str(aRight) ' + ' num2str(bRight) '  * (x^c)']);
fitdescriptionLeft = fittype([num2str(aLeft) ' + ' num2str(bLeft) '  * (x^c)']);
disp([fitdescriptionRight]);
disp([fitdescriptionLeft]);


%[foRight,gof] = fit(intendedLuminances', measuredLuminancesRightBlue', 'poly3');
%[foLeft, gof] = fit(intendedLuminances', measuredLuminancesLeftRed', 'poly3');
%x = linspace(0,1,256);
%measurmentSimulationRight = foRight.p1*x.^3 + foRight.p2*x.^2 + foRight.p3.*x + foRight.p4;
%measurmentSimulationLeft = foLeft.p1*x.^3   + foLeft.p2*x.^2  + foLeft.p3.*x  + foLeft.p4;


[foRight,gof,output] = fit(intendedLuminances', measuredLuminancesRightBlue', fitdescriptionRight);
[foLeft, gof,output] = fit(intendedLuminances', measuredLuminancesLeftRed', fitdescriptionLeft);
measuredGammaRight = foRight.c;
measuredGammaLeft = foLeft.c;
measurmentSimulationRight = aRight + bRight * linspace(0,1,256) .^measuredGammaRight;
measurmentSimulationLeft = aRight + bLeft * linspace(0,1,256) .^measuredGammaLeft;


figure(2182), clf, hold on, 
plot(linspace(0,1,256).*256, measurmentSimulationRight, '--c','LineWidth', 4);
plot(linspace(0,1,256).*256, measurmentSimulationLeft, '.-b', 'LineWidth', 4);
plot(intendedLuminances*256, measuredLuminancesRightBlue, 'm');
plot(intendedLuminances*256, measuredLuminancesLeftRed, 'r');

mycc = iccread(oldICCFilePath);

oldLookupTableRight = []; oldLookupTableLeft = [];
for i = 1:size(mycc.PrivateTags,1)
    if(strcmp(mycc.PrivateTags{i,1}, 'vcgt'))
        oldLookupTableRight = eight2sixteen(mycc.PrivateTags{i,2}(1043:1554)); %mycc.PrivateTags{i,2}(1043:2:1554);
        oldLookupTableLeft  = eight2sixteen(mycc.PrivateTags{i,2}(19:530)); %mycc.PrivateTags{i,2}(19:2:530);
        oldLookupTableRight = oldLookupTableRight ./ 256;
        oldLookupTableLeft = oldLookupTableLeft ./ 256;
        ptd = mycc.PrivateTags{i,2};
        vgctIndex = i;
        break
    end
end
if(isempty(oldLookupTableRight) | isempty(oldLookupTableLeft))
    disp('Cant read the lookuptable from ICC file');
    return;
end

measurmentSimulationScaledRight  = measurmentSimulationRight ./ max(measurmentSimulationRight) .* max(targetKernel); 
measurmentSimulationScaledLeft  = measurmentSimulationLeft   ./ max(measurmentSimulationLeft)  .* max(targetKernel); 
%yRight = oldLookupTableRight .* targetKernel ./ measurmentSimulationScaledRight;
%yLeft = oldLookupTableLeft   .* targetKernel ./ measurmentSimulationScaledLeft;

yRight = oldLookupTableRight .* (measurmentSimulationScaledRight ./ targetKernel);
yLeft = oldLookupTableLeft   .* (measurmentSimulationScaledLeft ./ targetKernel);

% for ii = 1:256
%     [min_difference, array_position] = min(abs(oldLookupTableRight - measurmentSimulationScaledRight(ii)));
%     yRight(ii) = measurmentSimulationScaledRight(array_position);
%     [min_difference, array_position] = min(abs(oldLookupTableLeft - measurmentSimulationScaledLeft(ii)));
%     yLeft(ii) = measurmentSimulationScaledLeft(array_position);
% end



figure(2919), clf, hold on
plot(double(oldLookupTableRight), 'b');
plot(double(oldLookupTableLeft), 'r');
plot(targetKernel, 'k');
plot(measurmentSimulationScaledRight, '--g', 'LineWidth', 3);
plot(measurmentSimulationScaledLeft, '.-g', 'LineWidth', 3);
plot(yRight, '*-m');
plot(yLeft, '.-m');

n16yRight = double(yRight);
n16yRight = (n16yRight - min(n16yRight)) ./ range(n16yRight);
n16yRight = n16yRight * 2^16;
n16yRight = uint16(n16yRight);
n16yRight = sixteen2eight(n16yRight);

n16yLeft = double(yLeft);
n16yLeft = (n16yLeft - min(n16yLeft)) ./ range(n16yLeft);
n16yLeft = n16yLeft * 2^16;
n16yLeft = uint16(n16yLeft);
n16yLeft = sixteen2eight(n16yLeft);



% 1:18 are the strings don't change them

% 19:530 are the fisrt gun shoud be RED
%ptd(19:2:530) = y; %high bits
%ptd(20:2:530) = 0; %low bits
ptd(19:530) = n16yLeft; 

% 531:1042 are the second gun shoud be GREEB
%ptd(531:2:1042) = y; %high bits
%ptd(532:2:1042) = 0; %low bits
ptd(531:1042) = mean([n16yLeft; n16yRight]);

% 1043:1554 are the third gun shoud be BLUE
%ptd(1043:2:1554) = y; %high bits
%ptd(1044:2:1554) = 0; %low bits
ptd(1043:1554) = n16yRight; 

mycc.PrivateTags{vgctIndex,2} = ptd;

% the name that appears in the preferences window in the list of profiles
%mycc.Description.String = ['ColorLCDGamma' num2str(measuredGamma)  '.', num2str(now)];
newName = ['Ali' char('A' + randi(24,10,1))'];
mycc.Description.String = newName;

rc = 1;
mycc.PrivateTags{1,2}(37) = uint8(0);
mycc.PrivateTags{1,2}(38) = uint8(0);
for ri = 1:2:35
    if(ri<=length(newName))
        mycc.PrivateTags{1,2}(39+ri) = uint8(newName(rc));
    else
        mycc.PrivateTags{1,2}(39+ri) = uint8(0);
    end
    rc = rc + 1;
end

iccwrite(mycc, [outputPath newName '.icc']);
results = ['success ' newName];

end

