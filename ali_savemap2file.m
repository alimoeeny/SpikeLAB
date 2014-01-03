function err = ali_savemap2file(m, filename)
% gets a map (containers.Map) and a filename (full path to the file) and
% saves a string representation of the map into the text file. this is
% primarily mean to to be used for binoc (stim) files, and puts each
% key-value pair in one line with a = between them without any spaces.
    mapstr = '';
    for k = keys(m), 
        mapstr = [mapstr, k{1}, '=', num2str(m(k{1})), '\n'];
    end
    fid = fopen(filename, 'w');
    fprintf(fid, mapstr);
    fclose(fid);
    err = '';
    