    function pos = PlaceUi(a, b ,str)        space = 0.01;        pos = get(b,'Position');        x = get(a,'position');        if strcmp(str,'left');            pos(1) = x(1)+x(3);            pos(2) = x(2);        elseif strcmp(str,'down');            pos(1) = x(1);            pos(2) = x(2)-x(4);        elseif strcmp(str,'up');            pos(1) = x(1);            pos(2) = x(2)+x(4);        end        set(b,'position', pos);        