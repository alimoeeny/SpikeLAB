winsize = 30;
threshold = 2;
swtchcounter = 1;
swatch = [];
lastdetection = 0;
for i = winsize: length(convspaceN) - winsize
    if(i>lastdetection + winsize)
        if (convspaceN(i)>threshold)
            pk = max(convspace(i:i+winsize));
            pkpoint = find(convspace(i:i+winsize)==pk, 1, 'first');
        
            swatch(swtchcounter,:) = V12(i + pkpoint - 20:i + pkpoint + 20);
            swtchcounter = swtchcounter + 1;
            lastdetection = i;
            disp(i);
        end
    end
end