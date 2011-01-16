function clickscatter(x, y, varargin)
% x , y, Color, size, filenames

mycolors = 'rgbcmykrgbcmykrgbcmykrgbcmykrgbcmyk';

srk = 1;
if (length(varargin) > 0)
    srk = varargin{1};
    if length(srk)<2
        srk = ones(length(x),1) + srk;
    end
else
    srk = (ones(length(x),1) + srk)';
end

srke = 1;
if (length(varargin) > 3)
    srke = varargin{4};
    if length(srke)<2
        srke = ones(length(x),1) + srke;
    end
else
    srke = (ones(length(x),1) + srke)';
end


sz = 7;
if (length(varargin) > 1)
    if ~isempty(varargin{2})
        if (length(varargin{2})>1)
            sz = varargin{2};
        else
            sz = ones(length(x),1) + varargin{2};
        end
    else
        sz = (ones(length(x),1) + sz)';
    end
end

if (length(varargin) > 2)
    if(~isempty(varargin{3}))
        if(length(varargin{3})>2)
            fileNames = varargin{3};
        else
            for i = 1:length(x), 
                fileNames{i} = num2str(i);
            end
        end
    else 
        for i = 1:length(x), 
            fileNames{i} = num2str(i);
        end
    end
end



hold on,
for i=1:length(x)
        plot(x(i), y(i),'-.or','MarkerEdgeColor', mycolors(srk(i)),...
                'MarkerFaceColor',mycolors(srk(i)), ...  %[.49 1 .63],...
                'MarkerEdgeColor',mycolors(srke(i)), ... 
                'MarkerSize',sz(i), ...
                'buttondownfcn',{@aliclicked, 'plotpsth', gcf, i, fileNames{i}});

end
refline(0,0);
reflinexy(0,100, 'LineStyle', '-');
hold off,
