function AddCellMenu(DATA)    if DATA.interactive < 0        return;    end    hm = findobj(DATA.toplevel,'Tag','ProbeMenu');    if ~isfield(DATA,'CellList')        return;    end    sm = findobj(hm, 'Tag','CellSwitchMenu');    delete(sm);    sm =  uimenu(hm,'Label','&Cell','Tag','CellSwitchMenu');    cellids = unique(DATA.CellList(:));    cellids = cellids(cellids > 0);    cellids = reshape(cellids,1,length(cellids));    pchars = ['1':'9' '0' 'a':'z'];    for j = cellids        if j < 10            h = uimenu(sm,'Label',['&' num2str(j)] ,'Callback',{@AllV.ChangeCell, j});        elseif j <= length(pchars)            uimenu(sm,'Label',[num2str(j) ' (&'  pchars(j) ')'] ,'Callback',{@AllV.ChangeCell, j});        else            uimenu(sm,'Label', num2str(j) ,'Callback',{@AllV.ChangeCell, j});        end    end        