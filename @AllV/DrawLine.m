function h= DrawLine(E,varargin)x = [E.pos(1) E.pos(3)];y = [E.pos(2) E.pos(4)];%fprintf('%.3f,%.3f -> %.3f,%.3f\n',x(1),y(1),x(2),y(2));if isfield(E,'h') & ishandle(E.h) & get(E.h,'parent') == E.axis;     set(E.h,'Xdata',x,'Ydata',y);    h = E.h;else    hold on;    h = plot(real(x),real(y),varargin{:});    hold off;end