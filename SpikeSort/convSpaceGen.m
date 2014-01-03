winsize = 30;
convspace = [];
parfor i = winsize:length(V12n)-winsize
    spliwin = V12n(i:i+winsize);
    convscapc(i) = sum(conv(spliwin, V12n, 'same'));
end

