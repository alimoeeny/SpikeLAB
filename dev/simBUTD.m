



TotalAverageFiringRate = 10; 
ButtomUpEffectSize = 3;      %firing rate change 
TopDownEffectSize = 3;       %firing rate change 


TotalTrialsCount = 1000;
TrialsWithOnlyButtomUpEffect = 500;
TrialsWithOnlyTopDownEffect = 500;

TrialCountsTableForBUTrials = [TrialsWithOnlyButtomUpEffect /2 TrialsWithOnlyButtomUpEffect /2 ; TrialsWithOnlyButtomUpEffect /2 TrialsWithOnlyButtomUpEffect /2 ];

TrialCountsTableForTDTrials = [TrialsWithOnlyTopDownEffect 0 ; TrialsWithOnlyTopDownEffect 0];



BU_PrefStimPrefChoiceTrials = TotalAverageFiringRate + ButtomUpEffectSize * randn(TrialCountsTableForBUTrials(1,1),1);

BU_PrefStimNullChoiseTrials = TotalAverageFiringRate - ButtomUpEffectSize * randn(TrialCountsTableForBUTrials(1,1),1);



%% descriptive

% tables = [choice+ stimulus+  choice- stimulus+; choice+ stimulus- choice- stimulus-];


%% Buttom up model

b = ButtomUpEffectSize;
trials_count_table = [n/4 n/4 ; n/4 n/4];
firing_rates_table = [b -b ; b , -b];



%% top down model

t = TopDownEffectSize;
trials_count_table = [n/2 0 ; 0 n/2];
firing_rates_table = [t NaN ; NaN , -t];

%% mixed model

b = ButtomUpEffectSize;
t = TopDownEffectSize;
trials_count_table = [n/4+m n/4-m ; n/4-m n/4+m];
firing_rates_table = [(m*t+b*(n/4))/(n/4+m) -b;  +b (-m*t-b*(n/4))/(n/4+m)];



