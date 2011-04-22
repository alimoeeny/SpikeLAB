function r = myParSave(filename, Expt)
    r = -1;
    try
        save(filename, 'Expt')
        r = 0;
    catch exception
        disp(exception.message);
        r = exception.identifier;
    end