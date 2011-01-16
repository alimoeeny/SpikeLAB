function [m, n, c] = NeurClus(NN)
if NN{1}(1) == 'd'
    NNN = str2num(NN{1}(4:end));
    m = 'dae';
else
    NNN = str2num(NN{1}(3:end));
    m = 'icarus';
end
if (NN{1}(1) == 'd' & str2num(NN{1}(4))==0)
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
