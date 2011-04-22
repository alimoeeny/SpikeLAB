function d = ali2Ddistance(a)
d = 0;
if size(a,1)>2 && size(a,2)==2
    for i = 2:size(a,1)
        d = d + sqrt((a(i,1)-a(i-1,1))^2 + (a(i,2)-a(i-1,2))^2); 
    end
else if size(a',1)>2 && size(a',2)==2
        d = ali2Ddistance(a');
    else
        d = -inf;
        disp('cant do the 2D distance input not in a good shape');
    end
end