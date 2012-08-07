function [pp, pn] = TailFits(prefPSTH, nullPSTH, selectRange)
% pp and pn could be fed to poly(p) to get the polynomial (here just a
% linear) function back
    
    [pp, sp, mup] = polyfit(selectRange, prefPSTH(selectRange),1);
    
    [pn, sn, mun] = polyfit(selectRange, nullPSTH(selectRange),1);
    
