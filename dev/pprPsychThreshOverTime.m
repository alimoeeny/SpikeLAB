clear, clc


flist = textread('pprpsychlist_all.txt', '%s');

Grand = {};
Grand{1} = []; % MIXED
Grand{2} = []; % -45
Grand{3} = []; % 0
Grand{4} = []; % 45
Grand{5} = []; % 90
Grand{6} = []; % MISSED

for i = 1:length(flist)
  disp([num2str(i) ' - ' flist{i}]);
  Expts = PsychMon(flist{i}, 'getexpts');
  for j = 1:length(Expts)
      Expt = Expts{j};
      if((length(Expt.Trials)>=60) && (length(unique([Expt.Trials.dx]))>4))
          s = PsychSlop(Expt, '', 'DID');
          if(~isempty(s))
              or = ExperimentProperties(Expt, 'or');
              try
                  if (length(unique([Expt.Trials.or]))>1)
                      or = -1;
                  end
              catch
              end
              switch or 
                  case -1
                      Grand{1}(end+1) = abs(s.fit(2));
                  case -45
                      Grand{2}(end+1) = abs(s.fit(2));
                  case 135
                      Grand{2}(end+1) = abs(s.fit(2));
                  case 0 
                      Grand{3}(end+1) = abs(s.fit(2));
                  case 45
                      Grand{4}(end+1) = abs(s.fit(2));
                  case 90
                      Grand{5}(end+1) = abs(s.fit(2));
                  otherwise
                      disp(['%$$%^#^$#*$%@(*#^#&^@&^#($!@        -> ' num2str(or)]);
                      Grand{6}(end+1) = abs(s.fit(2));
              end
          end
      end
  end
end

%%
figure(112), clf, hold on,
plot(Grand{1}, 'm')
plot(Grand{2}, 'b')
plot(Grand{3}, 'r')
plot(Grand{4}, 'g')
plot(Grand{5}, 'k')
plot(Grand{6}, 'c')
refline(0)


