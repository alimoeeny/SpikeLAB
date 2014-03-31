function hxy = reflinexy(x,y,varargin)
% REFLINEXY - Plot vertical and horizontal reference lines
%   REFLINEXY(X,Y) plots vertical and horizontal reference lines. The
%   vertical set of lines are drawn from the bottom of the axes upto each
%   point Y at the positions X. The horizontal set of lines are drawn from
%   the left of the axes to each point X at the positions Y.
%
%   REFLINEXY(X,Y,XOFFSET) uses the value XOFFSET as the vertical offset
%   for the vertical lines. REFLINE(X,Y,XOFFSET,YOFFSET) uses the value
%   YOFFSET as the horizontal offset for the horizontal lines. If the
%   offsets are omitted, or if they are empty, the defaults are used,
%   namely the bottom and/or the left side of the axis, respectively.
%   The offsets can be specified for each point separately.
%
%   The lines are plotted as a two graphics object. H = XYREFLINE(..) returns
%   the two graphics handles of the vertical and horizontal line objects. 
%
%   REFLINEXY(..., 'Prop1','Val1','Prop2','Val2', ...) uses the properties
%   and values specified for color, linestyle, etc. Execute GET(H), where H is
%   a line handle, to see a list of line object properties and their current values.
%   Execute SET(H) to see a list of line object properties and legal property values.
%   REFLINEXY uses the current axes, if any. Lines outside the plot area
%   are plotted but not shown.
%
%   Examples
%     % some data
%       N = 1000 ; i=0:N ;
%       x = (i/N) * 4 * pi ; y = exp(i/N) .* sin(x) ;
%       plot(x,y,'b-') ;
%     % add specific reference lines
%       i0 = [100 400 600] ;
%       reflinexy(x(i0),y(i0),'color','r') ;
%     % other reference lines with a specific offset
%       h = reflinexy(1:14,linspace(-2,1.5,14),[],x(end)-2) ;
%       set(h,'color','g','linestyle','-','linewidth',3) ;
%
%   See also STEM, PLOT, REFLINE, GRID, AXES
%   and GRIDXY (on the FEX)

% for Matlab R13
% version 1.2 (apr 2007)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0 (mar 2007) - created
% 1.1 (apr 2007) - included error check for line properties. Changed name,
%      and renamed XYRELINE.M to GRIDXY.M
% 1.2 (apr 2007) - added warning message idenitfier; changed first argument
%      to uistack to be a column vector (apparently something has changed
%      in uistack since R13)

error(nargchk(2,Inf,nargin)) ;

% check the arguments
if ~isnumeric(x) || ~isnumeric(y),
    error('Numeric arguments expected') ;
end

N = numel(x) ;

if numel(y) ~= N,
    error('X and Y should have the same number of elements') ;
end

% assume that the reference are to be draw ...
y0 = [] ; % ... from the bottom
x0 = [] ; % ... and from the left of the axis

if nargin==2,
    va = [] ;
else
    va = varargin ;
    if isnumeric(va{1})        
        % first optional argument specifies y offset        
        y0 = va{1} ;
        if ~isempty(y0),
            if numel(y0) == 1 ,
                y0 = repmat(y0,numel(x),1) ;
            elseif numel(y0) ~= N,
                error('All 4 matrices (x,y, offsetY, offsetX) should have the same number of elements.') ;
            end
        end
        va = va(2:end) ;
    elseif ~ischar(va{1})
        error('Invalid third argument') ;
    end
    if isnumeric(va{1})        
        % second optional argument specifies X offset        
        x0 = va{1} ;
        if ~isempty(x0),
            if numel(x0) == 1 ,
                x0 = repmat(x0,numel(x),1) ;
            elseif numel(x0) ~= N,
                error('All 4 matrices (x,y, offsetY, offsetX) should have the same number of elements.') ;
            end
            va = va(2:end) ;
        end
    elseif ~ischar(va{1})
        error('Invalid fourth argument') ;
    end
    if mod(size(va),2) == 1,
        error('Property-Value have to be pairs') ;
    end
end

% get the axes to plot in
hca=get(get(0,'currentfigure'),'currentaxes');
if isempty(hca),
    warning([mfilename ':NoAxes'],'No current axes found') ;
    return ;
end

% get the current limits of the axis
% also used for limit restoration later on
xlim = get(hca,'xlim') ;
ylim = get(hca,'ylim') ;

if isempty(y0), 
    % default: use bottom as the vertical offset
    y0 = repmat(ylim(1),N,1) ;
end
if isempty(x0),
    % default: use left as the horizontal offset
    x0 = repmat(xlim(1),N,1) ;
end

if ~isempty(x),   
    % setup data for the vertical lines
    xx1 = repmat(x(:).',3,1) ;
    yy1 = [y0(:)' ; y(:).' ; repmat(nan,1,numel(y))] ;
    % setup data for horizontal lines
    xx2 = [x0(:)' ; x(:).' ; repmat(nan,1,numel(y))] ;
    yy2 = repmat(y(:).',3,1) ;
    
    % add the line to the current axes
    np = get(hca,'nextplot') ;
    set(hca,'nextplot','add') ;
    hxy(1) = line('xdata',xx1(:),'ydata',yy1(:),'linestyle','-','color','k') ; % draw vertical lines
    hxy(2) = line('xdata',xx2(:),'ydata',yy2(:),'linestyle','-','color','k') ; % draw horizontal lines         
    
    uistack(hxy(:),'bottom') ; % push lines to the bottom of the graph
    
    set(hca,'nextplot',np) ;    % reset the nextplot state
    set(hca,'ylim',ylim,'xlim',xlim) ; % reset the limits
    
    if ~isempty(va),
        try
            set(hxy,va{:}) ; % set line properties        
        catch
            % invalid arguments, modify error message
            delete(hxy) ;
            msg = lasterror ;
            error(msg.message(21:end)) ;    
        end
    end
else
    hxy = [] ;
end

% if requested return handles
if ~nargout,     
    clear hxy ;
end









