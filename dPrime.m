function [dp] = dPrime(a,b, varargin)

if ~isempty(varargin)
    if strcmpi(varargin, 'r')
        temp = a;
        a = b;
        b = temp;
    end
end

if isempty(a) | isempty(b) 
    dp = -123456;
    disp('You''ve got to be kidding me, empty input for d prime!');
else

    aa = mean(a);
    bb = mean(b);
    av = var(a);
    bv = var(b);

    dp = (aa - bb) / sqrt((av + bv) /2);
end