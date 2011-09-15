% CP simulation

showplot = 0;

rcounter = 0;
for repetition = 1:10
    disp(repetition);
    
TrialsCount = 100;
% ->> choice1_samples = 0.5;
for choice1_samples = 0.0:0.05:1.0
    
choice2_samples = 1 - choice1_samples;

ScaleingFactor = 50;

% ->> knob = 0.2;
for knob = 0.05:0.05:0.5

choice1_1offset = knob * 1;
choice2_1offset = knob * -1;
choice1_2offset = knob * 1;
choice2_2offset = knob * -1;
choice_2_1shift = knob * 1.5;
choice_2_2shift = knob * -1.5;




%% Sampling 


choice1_1 = ScaleingFactor * (randn(TrialsCount,1) + choice1_1offset);
choice2_1 = ScaleingFactor * (randn(TrialsCount,1) + choice2_1offset);
choice1_2_1 = ScaleingFactor * (randn(TrialsCount * choice1_samples,1) + choice1_2offset + choice_2_1shift);
choice2_2_1 = ScaleingFactor * (randn(TrialsCount * choice2_samples,1) + choice2_2offset + choice_2_1shift);
choice1_2_2 = ScaleingFactor * (randn(TrialsCount * choice1_samples,1) + choice1_2offset + choice_2_2shift);
choice2_2_2 = ScaleingFactor * (randn(TrialsCount * choice2_samples,1) + choice2_2offset + choice_2_2shift);

if showplot
    figure(1111), clf, hold on,
    hist(choice1_1); 
    hist(choice2_1); 
    hist(choice1_2_1); 
    hist(choice2_2_1); 
    hist(choice1_2_2); 
    hist(choice2_2_2); 
end


%% CP at zero

CPatZero = ROCAUC(choice1_1, choice2_1);



%% CP with zscore

a = zscore([choice1_2_1; choice2_2_1]);
b = zscore([choice1_2_2; choice2_2_2]);

aa = [a(1:length(choice1_2_1)); b(1:length(choice1_2_2))];
bb = [a(length(choice1_2_1)+1 : length(choice1_2_1)+length(choice2_2_1)); b(length(choice1_2_2)+1 : length(choice1_2_2) + length(choice2_2_2))];

CPzscore = ROCAUC(aa, bb);


%% CP Weighted

CP_2_1 = ROCAUC(choice1_2_1, choice2_2_1);
CP_2_2 = ROCAUC(choice1_2_2, choice2_2_2);

CPW = 0.5 * ( CP_2_1 + CP_2_2 );


rcounter = rcounter + 1;
results(rcounter,:) = [choice1_samples , CPatZero, CPzscore, CPW, knob];

end % for knob = 0.05:0.05:0.5

end %  choice1_samples = 0.0:0.05:1.0

end % for repeatation = 1:10

%% Graphics

figure(2345), clf, hold on
scatter(results(:,1), results(:,2), 'r', 'filled'); refline(1)
scatter(results(:,1), results(:,3), 'b', 'filled'); refline(1)
scatter(results(:,1), results(:,4), 'g', 'filled'); refline(1)

figure(4567), clf, hold on
scatter(results(:,3), results(:,4), 'k', 'filled'); refline(1)
