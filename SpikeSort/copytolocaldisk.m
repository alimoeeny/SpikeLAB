tmp = load('../AllFullVFiles.mat', 'AllFullVFiles');
AllFullVFiles  = tmp.AllFullVFiles;
clear tmp;

for i = 1:length(AllFullVFiles)
   s = AllFullVFiles{i};
   a = strfind(s, '/');
   if(strcmpi('jbe', s(a(4)+1:a(5)-1))==0)
    s2 = ['/data' s(a(4):end)];
    if(exist(s, 'file')>0)
    if(exist(s2, 'file')==0)
        if(exist(['/data' s(a(4):a(6))], 'dir')==0)
            mkdir(['/data' s(a(4):a(6))]); 
            disp(['============    /data' s(a(4):a(6))  '--------> ']);
        end
        disp(s2);
        try
            copyfile(s, s2);
        catch err
            disp(err);
        end
    end
    else
        disp(['@#$#$#$@!! ==>   ' s '<<<<<<']);
    end
   end
end