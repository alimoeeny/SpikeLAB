function DATA = SetPCs(DATA, replot)    AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');    [C, DATA.Evec, pcs, dip, chspk, errs, DATA.pcfit] = AllV.CalcPCs(DATA, AllVoltages,DATA.nprobepc);    DATA.errs = {DATA.errs{:} errs{:}};    DATA.pcs(DATA.uid,1:size(pcs,2)) = pcs;    if replot    DATA = AllV.ReplotPCs(DATA,[]);    end