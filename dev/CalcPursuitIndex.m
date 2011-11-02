function [pursuiti, rates] = CalcPursuitIndex(Expt)


argon = {};
         Expt = SetPdir(Expt);
         Expt.Stimvals.e2 = 'pvel';
         if Expt.Stimvals.jx == 0
             dur = median([Expt.Trials.dur]);
             argon = {argon{:} 'latency'  500  'duration' dur};
         end
             
         thedata = PlotExpt(Expt,'nmin',5,argon{:});
         if isfield(thedata,'x')
%calculate prefered disp
         pid = find(thedata.x(:,1) > 0);
         nid = find(thedata.x(:,1) < 0);
         zid = find(thedata.x(:,1) == 0);
         if isempty(pid) || isempty(nid)
             dprime = NaN;
             dprimem = NaN;
         else
         for k = 1:length(nid)
         for j = 1:size(thedata.y,2)
             a = nid(k);
             b = pid(end-k+1);
             %dprime measure 
             sd = sqrt((thedata.sd(a,j).^2+thedata.sd(a,j).^2)./2);
             dprime(j,k) = (thedata.means(a,j) - thedata.means(b,j))./sd;
             dprimem(j,k) = dprime(j,k).*(sqrt(thedata.n(a,j)) + sqrt(thedata.n(b,j)));
         end
         end
         end
%dprime measures effect of disparities used on rate
%calculate effect of pursuit, for each disp
%dps is abs pursuit speed.
         dps = unique(abs(thedata.y(find(thedata.n > 0))));
         if isempty(pid) || isempty(nid) || isempty(zid)
             dprimes = NaN;
             rates = nan;
         else
         for k = 1:length(nid)
         for j = 1:length(dps)
             a = nid(k);
             b = pid(end-k+1);
             ip = find(thedata.y(1,:) == dps(j));
             in = find(thedata.y(1,:) == -dps(j));
             sd = sqrt((thedata.sd(a,ip).^2+thedata.sd(a,in).^2)./2);
             dprimes(j,1,k) = (thedata.means(a,ip) - thedata.means(a,in))./sd;
             sd = sqrt((thedata.sd(zid,ip).^2+thedata.sd(zid,in).^2)./2);
             dprimes(j,2,k) = (thedata.means(zid,ip) - thedata.means(zid,in))./sd;
             sd = sqrt((thedata.sd(b,ip).^2+thedata.sd(b,in).^2)./2);
             dprimes(j,3,k) = (thedata.means(b,ip) - thedata.means(b,in))./sd;
             rates(j,1,1) = thedata.means(a,ip);
             rates(j,1,2) = thedata.means(a,in);
             rates(j,2,1) = thedata.means(zid,ip);
             rates(j,2,2) = thedata.means(zid,in);
             rates(j,3,1) = thedata.means(b,ip);
             rates(j,3,2) = thedata.means(b,in);
         end
         end
         end
         pursuiti = diff(rates,[],3)./sum(rates,3);
         end
             
             pstr = [];
         if isempty(thedata) 
             title(sprintf('%s',splitpath(DATA.fstrings{id})));
         elseif isempty(thedata(1).title)
             title(sprintf('%s or%.0f',splitpath(thedata.name),GetEval(thedata.Data,'or')));
         else
             title(sprintf('%s or%.0f%s Fa %.0f',thedata(1).title,GetEval(thedata.Data,'or'),pstr,Expt.Stimvals.Pursuitdir));
         end
