
%% GENERATE A LINEAR ICC FILE
mycc = iccread('~/ColorLCD.icc');
pt = mycc.PrivateTags{4,2};
ptd = pt;
% 1:18 are the strings don't change them
% 19:530 are the fisrt gun shoud be RED
ptd(19:2:530) = [0:1:255]; %high bits
ptd(20:2:530) = 0; %low bits
% 531:1042 are the second gun shoud be GREEB
ptd(531:2:1042) = [0:1:255]; %high bits
ptd(532:2:1042) = 0; %low bits
% 1043:1554 are the third gun shoud be BLUE
ptd(1043:2:1554) = [0:1:255]; %high bits
ptd(1044:2:1554) = 0; %low bits


mycc.PrivateTags{4,2} = ptd;

iccwrite(mycc, '~/ColorLCDLinearHighabdLow.icc');



%% GENERATE AN ICC FILE WITH gamma =
gm = 1/2.2;
x = [1:1:256];
y = ((x/265).^gm);
y = ((y - min(y))./ (max(y)-min(y))) .* 256;

mycc = iccread('~/ColorLCD.icc');
pt = mycc.PrivateTags{4,2};
ptd = pt;
% 1:18 are the strings don't change them
% 19:530 are the fisrt gun shoud be RED
ptd(19:2:530) = y; %high bits
ptd(20:2:530) = 0; %low bits
% 531:1042 are the second gun shoud be GREEB
ptd(531:2:1042) = y; %high bits
ptd(532:2:1042) = 0; %low bits
% 1043:1554 are the third gun shoud be BLUE
ptd(1043:2:1554) = y; %high bits
ptd(1044:2:1554) = 0; %low bits


mycc.PrivateTags{4,2} = ptd;

mycc.Description.String = ['ColorLCDGamma' num2str(gm)];
iccwrite(mycc, ['~/ColorLCDGamma' num2str(gm) '.icc']);

