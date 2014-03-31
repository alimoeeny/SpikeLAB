function d = aliGazeEccentricity(a)
d = 0;
if size(a,1)>15 && size(a,2)==2
    ZeroX = mean(a(1:15,1));
    ZeroY = mean(a(1:15,2));
    for i = 16:size(a,1)
        d = d + sqrt((a(i,1)-ZeroX)^2 + (a(i,2)-ZeroY)^2); 
    end
else if size(a',1)>15 && size(a',2)==2
        d = aliGazeEccentricity(a');
    else
        d = -inf;
        disp('cant do the Gaze Eccentricity thing, input not in a good shape');
    end
end