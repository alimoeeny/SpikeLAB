function [deltafxy, Speeds] = dpiDeltaFXY(Expt)
        if isfield(Expt.Trials,'dfx')
            deltafxy = [Expt.Trials(:).dfx] - [Expt.Trials(:).fx];
            dfxy = [Expt.Trials(:).dfx];
            if ((Expt.Stimvals.or > 180) || (Expt.Stimvals.or < -180))
                deltafxy = -deltafxy;
                dfxy = - dfxy;
            end
        else
            deltafxy = [Expt.Trials(:).dfy] - [Expt.Trials(:).fy];
            dfxy = [Expt.Trials(:).dfy];
            if ((Expt.Stimvals.or < 90) && (Expt.Stimvals.or > -90))
                deltafxy = -deltafxy;
                dfxy = - dfxy;
            end
        end
        tempdelta = fix(deltafxy);
        tempd = fix(dfxy);
        for td = 1: length(tempdelta)
            if tempdelta(td) == 0
                tempdelta(td) = 0.5 * sign(deltafxy(td));
            end
        end
        deltafxy = tempdelta; % fix(deltafxy);  %  0.2 * fix(5*deltafxy); %0.1 * round(10*deltafxy);
        dfxy = tempd;
        Speeds = unique(abs(deltafxy));

