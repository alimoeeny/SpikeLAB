function x = GetValues(DATA, name)    x = [];    tid = find(strcmp(name,DATA.TemplateLabels));    if isempty(tid)        if sum(strncmp(name,{'ADC', 'Cen'},3))            AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');        end        if strncmp(name,'PC',2)           pc = sscanf(name,'PC%d');           x = DATA.pcs(:,pc);        elseif strncmp(name,'ADC',3)            if strfind(name,'dvdt')                AllVoltages = diff(AllVoltages,[],2);            end            if strncmp(name,'ADC1',4)                x = squeeze(AllVoltages(DATA.vpts(1,1),DATA.vpts(1,2),:));            elseif strncmp(name,'ADC2',4)                x =squeeze(AllVoltages(DATA.vpts(1,3),DATA.vpts(1,4),:));            end        elseif strcmp(name,'energy')            x = DATA.energy';        elseif strcmp(name,'spkvar')            x = DATA.spkvar';        elseif strcmp(name,'CentroidSplit')%Area under poisitive region post spike            t = find(DATA.spts ==0);            C = squeeze(AllVoltages(DATA.probe(1),t:end,:));            C(C < 0) = 0;            x = sum(C.^2)';        elseif strcmp(name,'CentroidSplit')%centroid of positiv points only                        C = squeeze(AllVoltages(DATA.probe(1),:,:));            C(C < 0) = 0;            x = C' * DATA.spts'./sum(C)';        elseif strcmp(name,'Centroid')            C = squeeze(AllVoltages(DATA.probe(1),:,:)).^2;            x = C'.^2 * DATA.spts'./sum(C.^2)';        end    else        x = DATA.TemplateScores(:,tid(1));    end            