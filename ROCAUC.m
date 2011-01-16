function [auc] = ROCAUC(a, b)

if isempty(a) || isempty(b) 
    auc = -1;
    return,
end
values = unique([a(:) ; b(:)]);
values(end+1) = max(values) + 1;
values(end+1) = min(values) - 1;
values = circshift(values,1);

for i = 1:length(values)
    TP = sum(a>values(i));
    TN = sum(b<=values(i));
    FP = sum(b>values(i));
    FN = sum(a<=values(i));
    xrocPoints(i) = 1 - (TN / (FP + TN));  %1- specificity
    yrocPoints(i) = TP / (TP + FN); % sensitivity
end

if(~issorted(xrocPoints))
    xrocPoints = flipud(xrocPoints');
    yrocPoints = flipud(yrocPoints');
end

    auc = 1-trapz(yrocPoints,xrocPoints);
if auc > 1
    disp('AUC > 1 !!! ');
end
end