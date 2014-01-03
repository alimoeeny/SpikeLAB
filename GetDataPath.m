function p = GetDataPath(~)
if(nargin>0)
    p = {'/Volumes/bgc/bgc/data/' '/Volumes/bgc6/smr/'};
else
    p = {'/bgc/data/' '/sd/bgc6/smr/'};
end