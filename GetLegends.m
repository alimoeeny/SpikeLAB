function leg = GetLegends(ft)
switch upper(ft)
    case 'TWO'
        leg = {'Pref bd Pref dx (all)', 'Pref bd null dx (all)', 'Pref bd ZERO dx', 'Pref bd ZERO dx (correct response)', 'null bd ZERO dx (correct response)', 'null bd ZERO dx', 'null bd Pref dx(flip)', 'Pref bd Pref dx', ' - ', 'pref bd null dx(flip)', 'null bd null dx'};
    case 'DID'
        leg = {'Preferred Id', 'Null Id', 'Pref Id Correct', 'Null Id Correct', 'Pref Id Wrong', 'Null Id Wrong', 'Pref dx', 'Null dx' }; %, 'Pref dx Top', 'Null dx Bottom'};
        %leg = {'Preferred Id', 'Null Id', 'Pref Id Correct', 'Null Id Correct', '-', 'Pref Id Wrong', '-', 'Null Id Wrong' }; %, 'Pref dx Top', 'Null dx Bottom'};
    case 'BDID'
        leg = {'Pref bd Pref Id', 'Pref bd null Id','Pref bd Pref Id P', 'Pref bd Pref Id and incorrect response', 'null bd Null Id C', 'null bd null Id incorrect response', 'Pref Id', 'null Id'};
    case 'RID'
        leg = {'Pref Or Pref dx','Pref Or null dx','Null Or Pref dx ','null Or null dx'};
    case 'DIDB'
        leg = {'Preferred Id', 'Null Id', 'Pref Id Correct', 'Null Id Correct', 'Pref Id Wrong', 'Null Id Wrong', 'Pref dx', 'Null dx' }; %, 'Pref dx Top', 'Null dx Bottom'};
    otherwise
        leg = {'Preferred Id', 'Null Id', 'Biased twoard pref', 'Biased toward Null'};
end
