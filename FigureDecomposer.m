function r = FigureDecomposer(h, filenamebase)

hlist = findobj(h, 'Type', 'line');
xl = get(h, 'xlim');
yl = get(h, 'ylim');

for i = 1: length(hlist)
    for j = 1: length(hlist)
        if (i~=j)
            set(hlist(j),'Visible','off');
        else
            set(hlist(j),'Visible','on');
        end
    end
    set(h, 'xlim', xl);
    set(h, 'ylim', yl);
    saveas(gcf, [filenamebase, num2str(i), '.pdf'], 'pdf');
end