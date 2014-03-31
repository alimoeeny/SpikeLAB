function [m, n, c] = NeurClus(NN)
switch NN{1}(1)
    case 'd'
        m = 'dae';
        if (NN{1}(length(m)+1)=='M')
            n = str2num(NN{1}(length(m)+2:end));
            c = [];
            return
        else
            NNN = str2num(NN{1}(4:end));
        end
        
    case 'i'
        m = 'icarus';
        NNN = str2num(NN{1}(3:end));
        
    case 't'
        m = 'test';
        NNN = str2num(NN{1}(4:end));
        
end

if (NN{1}(1) == 'd' & str2num(NN{1}(4))==0)
    if(NNN<100)
            n = NNN;
            c = '.c1.';
        else
            n = round(NNN/10);
            c = ['.c',num2str(mod(NNN,10)),'.'];
    end
elseif (NN{1}(1) == 'd' & str2num(NN{1}(4))==0)
    if(NNN<100)
            n = NNN;
            c = '.c1.';
        else
            n = round(NNN/10);
            c = ['.c',num2str(mod(NNN,10)),'.'];
    end    
else
    if(NNN<999)
            n = NNN;
            c = '.c1.';
        else
            n = round(NNN/10);
            c = ['.c',num2str(mod(NNN,10)),'.'];
    end
end
