function SaveCellList(DATA)    cellfile = [DATA.name '/CellList.mat'];    CellList = DATA.CellList;    CellDetails = DATA.CellDetails;    CellChanges = DATA.CellChanges;    save(cellfile, 'CellList','CellDetails','CellChanges');