function r = ali_randomizetrials(stimsdir)
	% find the index file and read it
    filename = ['/local/AliManStims/' stimsdir '/stimorder'];
    disp(['Randomizing: '  filename]);
    f = textread(filename);
    
    % randomize it
    g = f(randperm(length(f)));
    
    % save it back
    dlmwrite(filename, g, '\t');

end


