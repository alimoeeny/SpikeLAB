function DATA = CheckTemplates(DATA, C)if ~ismember(C.space(1),[3 4 6]) || isfield(DATA,'TemplateScores');    return;endDATA.TemplateLabels = AllV.TemplateLabels(DATA,0);Scores = AllV.CalcScores(DATA,C.MeanSpike);if size(Scores,1) > 1    DATA.TemplateScores(:,1)= Scores(2,1,:);    DATA.TemplateScores(:,8)= Scores(2,2,:);endDATA.TemplateScores(:,2)= sum(Scores(:,1,:));DATA.TemplateScores(:,3)= Scores(1,1,:);if size(Scores,1) > 2    DATA.TemplateScores(:,4)= Scores(3,1,:);endDATA.TemplateScores(:,10)= sum(Scores(:,2,:));DATA.TemplateScores(:,12)= 0;if DATA.trigdt == 4    DATA.TemplateScores(:,11) =DATA.rV;else    DATA.TemplateScores(:,11) = sum(TemplateScores(:,3,:));endDATA.tmpdips = AllV.CalculateTemplateDips(DATA);